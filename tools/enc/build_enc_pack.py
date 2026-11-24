#!/usr/bin/env python3
"""Build a Tiller Time vector chart pack (HENC) from free NOAA S-57 ENC cells.

NOAA ENCs are public domain: https://charts.noaa.gov/ENCs/ENCs.shtml
The pack keeps only what the app draws: depth areas, land, contours, coastline,
soundings, hazards, buoys, beacons, lights and place names.

Layout: b"HENC" + u32 header length + JSON header + raw DEFLATE body.
Body: a run of records (little endian):
  u8 kind, u8 band, u8 minZoom, u8 coveredBy, u32 attr, f32 v1, f32 v2,
  u8 labelLen, label (utf-8), varint parts, then per part: varint n, n x (zigzag-varint lon, lat)
  in 1e-6 deg, deltas from the previous point (the first from 0); kind 10 (soundings) adds
  n x zigzag-varint depth in decimetres after its one part.
Not for navigation.
"""
import argparse, json, math, os, re, struct, sys, urllib.request, zipfile, zlib, concurrent.futures as cf
from osgeo import ogr, gdal

gdal.UseExceptions()
gdal.SetConfigOption("OGR_S57_OPTIONS", "SPLIT_MULTIPOINT=ON,ADD_SOUNDG_DEPTH=ON,RETURN_PRIMITIVES=OFF,LNAM_REFS=OFF,RECODE_BY_DSSI=ON")

CATALOG = "https://charts.noaa.gov/ENCs/ENCProdCat_19115.xml"
CELL_URL = "https://charts.noaa.gov/ENCs/{}.zip"

K_COVER, K_DEPARE, K_DRGARE, K_LAND, K_LAKE, K_STRUCT, K_STRUCTLINE, K_BUILT, K_DEPCNT, K_COAST, \
    K_SOUND, K_HAZARD, K_HAZAREA, K_NAVAID, K_LIGHT, K_NAME = range(16)

NAVAID_OBJ = {"BOYLAT": 1, "BOYCAR": 2, "BOYISD": 3, "BOYSAW": 4, "BOYSPP": 5,
              "BCNLAT": 6, "BCNCAR": 7, "BCNISD": 8, "BCNSAW": 9, "BCNSPP": 10}
LITCHR = {1: "F", 2: "Fl", 3: "LFl", 4: "Q", 5: "VQ", 6: "UQ", 7: "Iso", 8: "Oc", 9: "IQ", 10: "IVQ",
          11: "IUQ", 12: "Mo", 13: "FFl", 14: "Fl+LFl", 15: "OcFl", 16: "FLFl", 25: "Q+LFl", 26: "VQ+LFl",
          27: "UQ+LFl", 28: "Al"}
COLOUR = {1: "W", 3: "R", 4: "G", 6: "Y", 5: "Bu", 11: "Or", 2: "Bk"}


def bbox_filter(cells, box):
    w, s, e, n = box
    return {k: v for k, v in cells.items() if v["e"] > w and v["w"] < e and v["n"] > s and v["s"] < n}


def read_catalog(path_or_url):
    data = open(path_or_url).read() if os.path.exists(path_or_url) else \
        urllib.request.urlopen(urllib.request.Request(path_or_url, headers={"User-Agent": "TillerTime"}), timeout=120).read().decode()
    out = {}
    for blk in data.split("<composedOf>")[1:]:
        t = re.search(r"<title>\s*<gco:CharacterString>([^<]+)", blk)
        pts = [tuple(map(float, p.split())) for p in re.findall(r"<gml:pos>([^<]+)</gml:pos>", blk)]
        if not t or not pts:
            continue
        a = [p[0] for p in pts]; b = [p[1] for p in pts]
        lat, lon = (b, a) if max(abs(x) for x in a) > 90 else (a, b)
        out[t.group(1)] = dict(w=min(lon), e=max(lon), s=min(lat), n=max(lat))
    return out


def fetch(name, workdir):
    path = os.path.join(workdir, name + ".000")
    if os.path.exists(path):
        return path
    z = os.path.join(workdir, name + ".zip")
    for attempt in range(3):
        try:
            urllib.request.urlretrieve(CELL_URL.format(name), z)
            break
        except Exception:
            if attempt == 2:
                raise
    with zipfile.ZipFile(z) as zf:
        for m in zf.namelist():
            if m.endswith(name + ".000"):
                open(path, "wb").write(zf.read(m))
    os.remove(z)
    return path


def scamin_zoom(f):
    """S-57 SCAMIN (1:N) to the first MapKit zoom showing it (3x iPhone, ~48 N)."""
    try:
        v = f.GetField("SCAMIN")
    except Exception:
        v = None
    if not v:
        return 0
    return max(0, min(22, math.ceil(math.log2(2.97e8 / float(v)))))


def field(f, name, default=None):
    i = f.GetFieldIndex(name)
    if i < 0 or not f.IsFieldSet(i):
        return default
    return f.GetField(i)


def first_int(v):
    if v is None:
        return 0
    if isinstance(v, list):
        return int(v[0]) if v else 0
    try:
        return int(str(v).split(",")[0])
    except ValueError:
        return 0


def ints(v):
    if v is None:
        return []
    if isinstance(v, list):
        return [int(x) for x in v if str(x).strip()]
    return [int(x) for x in str(v).split(",") if x.strip()]


TOLERANCE = {3: 1.5e-4, 4: 4e-5, 5: 1e-5}   # degrees: ~11 m, 3 m, 0.7 m


def simplified(geom, band):
    t = ogr.GT_Flatten(geom.GetGeometryType())
    if t in (ogr.wkbPoint, ogr.wkbMultiPoint):
        return geom
    g = geom.SimplifyPreserveTopology(TOLERANCE[band])
    return g if g is not None and not g.IsEmpty() else geom


def parts_of(geom, band=None):
    """Return list of coordinate lists (lon, lat) for any geometry."""
    if band is not None:
        geom = simplified(geom, band)
    out = []
    t = ogr.GT_Flatten(geom.GetGeometryType())
    if t in (ogr.wkbPolygon, ogr.wkbLineString, ogr.wkbLinearRing, ogr.wkbPoint):
        if t == ogr.wkbPolygon:
            for i in range(geom.GetGeometryCount()):
                out.append(geom.GetGeometryRef(i).GetPoints())
        else:
            out.append(geom.GetPoints())
    else:
        for i in range(geom.GetGeometryCount()):
            out += parts_of(geom.GetGeometryRef(i))
    return [[(p[0], p[1]) for p in pts] for pts in out if pts]


class Writer:
    def __init__(self):
        self.buf = bytearray()
        self.count = {}

    def var(self, v):
        while v >= 0x80:
            self.buf.append((v & 0x7F) | 0x80)
            v >>= 7
        self.buf.append(v)

    def zz(self, v):
        self.var((v << 1) ^ (v >> 63) if v < 0 else v << 1)

    def rec(self, kind, band, minz, attr, v1, v2, label, parts, depths=None, covered=0):
        label = (label or "").encode("utf-8")[:255]
        self.buf += struct.pack("<BBBBIffB", kind, band, minz, covered, attr & 0xFFFFFFFF,
                                float("nan") if v1 is None else v1, float("nan") if v2 is None else v2, len(label))
        self.buf += label
        self.var(len(parts))
        for pts in parts:
            self.var(len(pts))
            px = py = 0
            for x, y in pts:
                ix, iy = round(x * 1e6), round(y * 1e6)
                self.zz(ix - px); self.zz(iy - py)
                px, py = ix, iy
        if depths is not None:
            for d in depths:
                self.zz(round(d * 10))
        self.count[kind] = self.count.get(kind, 0) + 1


def light_text(f):
    ch = LITCHR.get(first_int(field(f, "LITCHR")), "")
    grp = field(f, "SIGGRP") or ""
    grp = grp.replace("()", "").strip()
    cols = "".join(COLOUR.get(c, "") for c in ints(field(f, "COLOUR")))
    per = field(f, "SIGPER")
    s = ch + (grp if grp and grp != "(1)" else "")
    s += (" " + cols) if cols else ""
    if per:
        s += " %gs" % float(per)
    return s.strip()


def short_name(name):
    if not name:
        return ""
    m = re.search(r'"?([0-9]+[A-Z]?|[A-Z]{1,2})"?$', name.strip())
    return m.group(1) if m else ""


def process(path, band, w, coverage):
    ds = ogr.Open(path)
    if ds is None:
        return

    def layer(name):
        return ds.GetLayerByName(name)

    def covered_by(geom):
        best = 0
        for b, g in coverage:
            if b > band and b > best and g.Contains(geom):
                best = b
        return best

    lyr = layer("M_COVR")
    if lyr:
        for f in lyr:
            if first_int(field(f, "CATCOV")) == 1 and f.GetGeometryRef():
                w.rec(K_COVER, band, 0, 0, None, None, "", parts_of(f.GetGeometryRef(), band))

    def areas(name, kind, v1f=None, v2f=None, polys_only=True):
        lyr = layer(name)
        if not lyr:
            return
        for f in lyr:
            g = f.GetGeometryRef()
            if g is None:
                continue
            gt = ogr.GT_Flatten(g.GetGeometryType())
            if polys_only and gt not in (ogr.wkbPolygon, ogr.wkbMultiPolygon):
                continue
            v1 = field(f, v1f) if v1f else None
            v2 = field(f, v2f) if v2f else None
            w.rec(kind, band, scamin_zoom(f), 0, v1, v2, "", parts_of(g, band))

    areas("DEPARE", K_DEPARE, "DRVAL1", "DRVAL2")
    areas("DRGARE", K_DRGARE, "DRVAL1", "DRVAL2")
    areas("LNDARE", K_LAND)
    areas("LAKARE", K_LAKE)
    areas("RIVERS", K_LAKE)
    areas("BUAARE", K_BUILT)
    for nm in ("SLCONS", "PONTON", "DRYDOC", "FLODOC", "HULKES", "BUISGL"):
        lyr = layer(nm)
        if not lyr:
            continue
        for f in lyr:
            g = f.GetGeometryRef()
            if g is None:
                continue
            gt = ogr.GT_Flatten(g.GetGeometryType())
            if gt in (ogr.wkbPolygon, ogr.wkbMultiPolygon):
                w.rec(K_STRUCT, band, scamin_zoom(f), 0, None, None, "", parts_of(g, band))
            elif gt in (ogr.wkbLineString, ogr.wkbMultiLineString):
                w.rec(K_STRUCTLINE, band, scamin_zoom(f), 0, None, None, "", parts_of(g, band))
    lyr = layer("BRIDGE")
    if lyr:
        for f in lyr:
            g = f.GetGeometryRef()
            if g is not None and ogr.GT_Flatten(g.GetGeometryType()) in (ogr.wkbLineString, ogr.wkbMultiLineString):
                w.rec(K_STRUCTLINE, band, scamin_zoom(f), 1, field(f, "VERCLR"), None, "", parts_of(g, band))

    lyr = layer("DEPCNT")
    if lyr:
        for f in lyr:
            if f.GetGeometryRef():
                w.rec(K_DEPCNT, band, scamin_zoom(f), 0, field(f, "VALDCO"), None, "", parts_of(f.GetGeometryRef(), band))
    lyr = layer("COALNE")
    if lyr:
        for f in lyr:
            if f.GetGeometryRef():
                w.rec(K_COAST, band, scamin_zoom(f), 0, None, None, "", parts_of(f.GetGeometryRef(), band))

    # Soundings: one record per (band, minZoom, coveredBy) group.
    lyr = layer("SOUNDG")
    if lyr:
        groups = {}
        for f in lyr:
            g = f.GetGeometryRef()
            if g is None:
                continue
            d = field(f, "DEPTH")
            if d is None:
                d = g.GetZ()
            key = (scamin_zoom(f), covered_by(g))
            groups.setdefault(key, ([], []))
            groups[key][0].append((g.GetX(), g.GetY()))
            groups[key][1].append(float(d))
        for (mz, cov), (pts, ds_) in groups.items():
            w.rec(K_SOUND, band, mz, 0, None, None, "", [pts], ds_, covered=cov)

    for nm, typ in (("UWTROC", 1), ("WRECKS", 2), ("OBSTRN", 3)):
        lyr = layer(nm)
        if not lyr:
            continue
        for f in lyr:
            g = f.GetGeometryRef()
            if g is None:
                continue
            attr = typ | (first_int(field(f, "WATLEV")) << 4) | (first_int(field(f, "CATWRK")) << 8) | \
                (first_int(field(f, "CATOBS")) << 12)
            val = field(f, "VALSOU")
            gt = ogr.GT_Flatten(g.GetGeometryType())
            if gt == ogr.wkbPoint:
                w.rec(K_HAZARD, band, scamin_zoom(f), attr, val, None, "", parts_of(g, band), covered=covered_by(g))
            elif gt in (ogr.wkbPolygon, ogr.wkbMultiPolygon, ogr.wkbLineString):
                w.rec(K_HAZAREA, band, scamin_zoom(f), attr, val, None, "", parts_of(g, band))

    for nm, code in NAVAID_OBJ.items():
        lyr = layer(nm)
        if not lyr:
            continue
        for f in lyr:
            g = f.GetGeometryRef()
            if g is None or ogr.GT_Flatten(g.GetGeometryType()) != ogr.wkbPoint:
                continue
            cat = first_int(field(f, "CATLAM")) or first_int(field(f, "CATCAM")) or first_int(field(f, "CATSPM"))
            shp = first_int(field(f, "BOYSHP")) or first_int(field(f, "BCNSHP"))
            cols = ints(field(f, "COLOUR"))[:3] + [0, 0, 0]
            attr = code | (min(cat, 15) << 4) | (min(shp, 15) << 8) | (cols[0] << 12) | (cols[1] << 16) | (cols[2] << 20) \
                | (min(first_int(field(f, "COLPAT")), 15) << 24)
            w.rec(K_NAVAID, band, scamin_zoom(f), attr, None, None, short_name(field(f, "OBJNAM")), parts_of(g, band),
                  covered=covered_by(g))

    lyr = layer("LIGHTS")
    if lyr:
        for f in lyr:
            g = f.GetGeometryRef()
            if g is None:
                continue
            col = first_int(field(f, "COLOUR"))
            attr = col | (min(first_int(field(f, "LITCHR")), 63) << 4)
            w.rec(K_LIGHT, band, scamin_zoom(f), attr, field(f, "SECTR1"), field(f, "SECTR2"), light_text(f),
                  parts_of(g, band), covered=covered_by(g))

    for nm in ("SEAARE", "LNDRGN"):
        lyr = layer(nm)
        if not lyr:
            continue
        for f in lyr:
            g = f.GetGeometryRef()
            name = field(f, "OBJNAM")
            if g is None or not name:
                continue
            p = g if ogr.GT_Flatten(g.GetGeometryType()) == ogr.wkbPoint else g.PointOnSurface()
            if p is None:
                continue
            w.rec(K_NAME, band, max(scamin_zoom(f), 10), 1 if nm == "SEAARE" else 2, None, None, name,
                  [[(p.GetX(), p.GetY())]], covered=covered_by(p))


def coverage_of(path, band):
    ds = ogr.Open(path)
    lyr = ds.GetLayerByName("M_COVR") if ds else None
    out = []
    if lyr:
        for f in lyr:
            if first_int(field(f, "CATCOV")) == 1 and f.GetGeometryRef():
                out.append((band, f.GetGeometryRef().Clone()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bbox", default="-125.3,47.0,-122.0,49.4", help="west,south,east,north")
    ap.add_argument("--bands", default="345")
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--work", default="enc_cells")
    ap.add_argument("--output", default="charts_salish.ttenc")
    a = ap.parse_args()
    box = tuple(map(float, a.bbox.split(",")))
    os.makedirs(a.work, exist_ok=True)
    cells = bbox_filter(read_catalog(a.catalog), box)
    names = sorted(k for k in cells if k[2] in a.bands)
    print("cells:", len(names), file=sys.stderr)
    with cf.ThreadPoolExecutor(8) as ex:
        paths = dict(zip(names, ex.map(lambda n: fetch(n, a.work), names)))
    coverage = []
    for n in names:
        coverage += coverage_of(paths[n], int(n[2]))
    w = Writer()
    for n in sorted(names, key=lambda k: (k[2], k)):
        band = int(n[2])
        near = [(b, g) for b, g in coverage if b > band]
        process(paths[n], band, w, near)
    body = zlib.compress(bytes(w.buf), 9)[2:-4]
    import datetime
    header = json.dumps({"magic": "HENC", "version": 1, "source": "NOAA ENC (public domain)",
                         "generated": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
                         "bbox": box, "cells": len(names), "records": sum(w.count.values()),
                         "counts": {str(k): v for k, v in sorted(w.count.items())},
                         "rawBytes": len(w.buf)}).encode()
    with open(a.output, "wb") as f:
        f.write(b"HENC" + struct.pack("<I", len(header)) + header + body)
    print("records", w.count, "raw %.1f MB, packed %.1f MB" % (len(w.buf) / 1e6, os.path.getsize(a.output) / 1e6),
          file=sys.stderr)


if __name__ == "__main__":
    main()
