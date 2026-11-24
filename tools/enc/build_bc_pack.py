#!/usr/bin/env python3
"""Build Tiller Time's free BC vector chart pack (HENC v2, band 6) from open data.

Sources (all free):
  * Depths: CHS NONNA-10 bathymetry (10 m grid), Open Government Licence - Canada,
    read from CHS's public WCS. Not for navigation.
  * Land: OpenStreetMap land polygons (osmdata.openstreetmap.de), ODbL.
  * Buoys, beacons, lights, rocks, wrecks, docks and place names: OpenStreetMap via Overpass, ODbL.

Writes the same record layout as build_enc_pack.py so the app draws both with one renderer.
The app clips this pack to the Canadian side of the border.

CHS's WCS refuses requests from GitHub's cloud runners, so the build has two stages:
  1. --depths-only (run on any ordinary machine): NONNA depths -> data/bc_depths.henc, committed
     to the repo; the survey data changes slowly, so rebuild it a few times a year.
  2. --depths-from bc_depths.henc (GitHub Actions, monthly): adds land, buoys, lights and names
     from OpenStreetMap and writes charts_bc.ttenc.
"""
import argparse, concurrent.futures as cf, datetime, json, math, os, re, struct, sys, tempfile, time, urllib.parse, \
    urllib.request, zipfile, zlib

import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
ogr.UseExceptions()

sys.path.insert(0, os.path.dirname(__file__))
from build_enc_pack import (Writer, merc, K_COVER, K_DEPARE, K_LAND, K_STRUCT, K_STRUCTLINE, K_DEPCNT, K_SOUND,
                            K_HAZARD, K_NAVAID, K_LIGHT, K_NAME, WORLD)

BAND = 6
R = 20037508.342789244
WCS = ("https://nonna-geoserver.data.chs-shc.ca/geoserver/wcs?service=WCS&version=1.0.0&request=GetCoverage"
       "&coverage=nonna:NONNA%2010%20Coverage&crs=EPSG:3857&response_crs=EPSG:3857"
       "&bbox={},{},{},{}&width={}&height={}&format=GeoTIFF")
LAND_URL = "https://osmdata.openstreetmap.de/download/land-polygons-split-4326.zip"
OVERPASS = ["https://overpass-api.de/api/interpreter", "https://overpass.kumi.systems/api/interpreter"]

PX = 10.0                 # metres (Web Mercator) per grid cell
CHUNK = 2048              # cells per chunk side (20 km)
AREA_LEVELS = [-20, -10, -5, -2, 0]                  # elevation (m, negative = depth), ascending
LINE_LEVELS = [-200, -100, -50, -30, -20, -10, -5, -2]
SOUND_CELLS = [(1600, 10), (800, 11), (400, 12), (200, 13), (100, 14), (70, 15)]   # cell size (m), first zoom
SIMPLIFY = 2.5            # metres


def m2w(x, y):
    """Web Mercator metres to pack world units."""
    return round((x + R) / (2 * R) * WORLD), round((R - y) / (2 * R) * WORLD)


def lonlat_to_m(lon, lat):
    return lon * R / 180, math.log(math.tan((90 + lat) * math.pi / 360)) * R / math.pi


def get(url, path=None, tries=4, timeout=300):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "TillerTime chart builder"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
            if path:
                open(path, "wb").write(data)
            return data
        except Exception as e:
            if k == tries - 1:
                raise
            time.sleep(5 * (k + 1))


class MWriter(Writer):
    """Writer that takes geometry already in pack world units."""

    def recw(self, kind, minz, attr, v1, v2, label, parts, depths=None):
        label = (label or "").encode("utf-8")[:255]
        self.buf += struct.pack("<BBBBIffB", kind, BAND, minz, 0, attr & 0xFFFFFFFF,
                                float("nan") if v1 is None else v1, float("nan") if v2 is None else v2, len(label))
        self.buf += label
        self.var(len(parts))
        for pts in parts:
            self.var(len(pts))
            px = py = 0
            for ix, iy in pts:
                self.zz(ix - px); self.zz(iy - py)
                px, py = ix, iy
        if depths is not None:
            for d in depths:
                self.zz(round(d * 10))
        self.count[kind] = self.count.get(kind, 0) + 1


def geom_parts_m(g):
    """Parts of an OGR geometry in metres -> world units, dropping repeated points."""
    out = []
    t = ogr.GT_Flatten(g.GetGeometryType())
    if t in (ogr.wkbPolygon,):
        for i in range(g.GetGeometryCount()):
            out += geom_parts_m(g.GetGeometryRef(i))
        return out
    if t in (ogr.wkbMultiPolygon, ogr.wkbMultiLineString, ogr.wkbGeometryCollection):
        for i in range(g.GetGeometryCount()):
            out += geom_parts_m(g.GetGeometryRef(i))
        return out
    pts = []
    for p in g.GetPoints() or []:
        w = m2w(p[0], p[1])
        if not pts or pts[-1] != w:
            pts.append(w)
    return [pts] if len(pts) > 1 else []


def simplify(g, preserve=True):
    try:
        r = g.SimplifyPreserveTopology(SIMPLIFY) if preserve else g.Simplify(SIMPLIFY)
        return r if r is not None and not r.IsEmpty() else g
    except RuntimeError:
        return g


def blur(a, mask):
    """3x3 mean of valid cells (light smoothing for tidy contours)."""
    v = np.where(mask, a, 0.0).astype(np.float64)
    m = mask.astype(np.float64)
    pv = np.pad(v, 1); pm = np.pad(m, 1)
    s = sum(pv[1 + dy:pv.shape[0] - 1 + dy, 1 + dx:pv.shape[1] - 1 + dx] for dy in (-1, 0, 1) for dx in (-1, 0, 1))
    n = sum(pm[1 + dy:pm.shape[0] - 1 + dy, 1 + dx:pm.shape[1] - 1 + dx] for dy in (-1, 0, 1) for dx in (-1, 0, 1))
    out = np.where(n > 0, s / np.maximum(n, 1), 0.0)
    return out


def contour(ds_mem, levels, polygons):
    drv = ogr.GetDriverByName("Memory")
    src = drv.CreateDataSource("c")
    srs = osr.SpatialReference(); srs.ImportFromEPSG(3857)
    lyr = src.CreateLayer("c", srs, ogr.wkbPolygon if polygons else ogr.wkbLineString)
    if polygons:
        lyr.CreateField(ogr.FieldDefn("emin", ogr.OFTReal)); lyr.CreateField(ogr.FieldDefn("emax", ogr.OFTReal))
        opts = ["FIXED_LEVELS=" + ",".join(str(l) for l in levels), "ELEV_FIELD_MIN=0", "ELEV_FIELD_MAX=1",
                "POLYGONIZE=YES", "NODATA=-99999"]
    else:
        lyr.CreateField(ogr.FieldDefn("elev", ogr.OFTReal))
        opts = ["FIXED_LEVELS=" + ",".join(str(l) for l in levels), "ELEV_FIELD=0", "NODATA=-99999"]
    gdal.ContourGenerateEx(ds_mem.GetRasterBand(1), lyr, options=opts)
    return src, lyr


LAND_LAYER = None
LAND_LOCK = None


def load_land_layer(box, zip_path):
    """OSM land in Web Mercator, in memory, for masking soundings that fall ashore."""
    import threading
    global LAND_LAYER, LAND_LOCK
    shp = "/vsizip/" + zip_path + "/land-polygons-split-4326/land_polygons.shp"
    src_ds = ogr.Open(shp)
    lyr = src_ds.GetLayer(0)
    lyr.SetSpatialFilterRect(*box)
    s4326 = osr.SpatialReference(); s4326.ImportFromEPSG(4326)
    s4326.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    s3857 = osr.SpatialReference(); s3857.ImportFromEPSG(3857)
    tr = osr.CoordinateTransformation(s4326, s3857)
    mem = ogr.GetDriverByName("Memory").CreateDataSource("land")
    out = mem.CreateLayer("land", s3857, ogr.wkbMultiPolygon)
    for f in lyr:
        g = f.GetGeometryRef().Clone()
        g.Transform(tr)
        nf = ogr.Feature(out.GetLayerDefn()); nf.SetGeometry(g); out.CreateFeature(nf)
    LAND_LAYER = (mem, out)
    LAND_LOCK = threading.Lock()


def land_mask(gt, shape):
    if LAND_LAYER is None:
        return np.zeros(shape, bool)
    r = gdal.GetDriverByName("MEM").Create("", shape[1], shape[0], 1, gdal.GDT_Byte)
    r.SetGeoTransform(gt)
    srs = osr.SpatialReference(); srs.ImportFromEPSG(3857); r.SetProjection(srs.ExportToWkt())
    with LAND_LOCK:
        lyr = LAND_LAYER[1]
        lyr.SetSpatialFilterRect(gt[0], gt[3] + shape[0] * gt[5], gt[0] + shape[1] * gt[1], gt[3])
        gdal.RasterizeLayer(r, [1], lyr, burn_values=[1], options=["ALL_TOUCHED=TRUE"])
        lyr.SetSpatialFilter(None)
    return r.GetRasterBand(1).ReadAsArray().astype(bool)


def process_chunk(args):
    x0, y1, work = args              # top-left corner in metres
    n = CHUNK + 1                    # one cell of overlap with the next chunk
    x1, y0 = x0 + n * PX, y1 - n * PX
    path = os.path.join(work, "n_%d_%d.tif" % (int(x0), int(y1)))
    if not os.path.exists(path):
        get(WCS.format(x0, y0, x1, y1, n, n), path)
    ds = gdal.Open(path)
    a = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    valid = np.isfinite(a) & (a < 1e30) & (a > -12000)
    if valid.sum() < 50:
        return None
    gt = ds.GetGeoTransform()
    # Close pinholes and survey-line gaps (up to ~30 m) so the tints don't speckle; land is
    # drawn from OpenStreetMap on top, so growing water under the shoreline is harmless.
    filled, fmask = a.copy(), valid.copy()
    for _ in range(3):
        b = blur(filled, fmask)
        grow = ~fmask & (blur(fmask.astype(float), np.ones_like(fmask)) > 0)
        filled[grow] = b[grow]
        fmask = fmask | grow
    sm = filled
    for _ in range(4):                        # ~15 m smoothing: tidy contours over blocky patches
        sm = blur(sm, fmask)
    sm[~fmask] = -99999
    mem = gdal.GetDriverByName("MEM").Create("", a.shape[1], a.shape[0], 1, gdal.GDT_Float32)
    mem.SetGeoTransform(gt)
    mem.GetRasterBand(1).WriteArray(sm.astype(np.float32))
    mem.GetRasterBand(1).SetNoDataValue(-99999)

    out = {"areas": [], "lines": [], "sound": []}
    src, lyr = contour(mem, AREA_LEVELS, True)
    for f in lyr:
        g = f.GetGeometryRef()
        if g is None:
            continue
        emin, emax = f.GetField(0), f.GetField(1)
        # The band's shallow limit as depth; anything above 0 m elevation dries.
        if emax > 1e6:
            d1 = -1.0                         # top interval (above 0 m): dries
        else:
            d1 = -emax                        # e.g. -5..-2 m elevation -> 2 m shallow limit; deepest -> 20
        g = simplify(g)
        if g.GetGeometryType() in (ogr.wkbMultiPolygon, ogr.wkbMultiPolygon25D):
            keep = ogr.Geometry(ogr.wkbMultiPolygon)
            for i in range(g.GetGeometryCount()):
                if g.GetGeometryRef(i).GetArea() >= 6000:     # isolated blocks under ~80 x 80 m are survey patchwork
                    keep.AddGeometry(g.GetGeometryRef(i))
            g = keep
        parts = geom_parts_m(g)
        if parts:
            out["areas"].append((d1, parts))
    src, lyr = contour(mem, LINE_LEVELS, False)
    for f in lyr:
        g = f.GetGeometryRef()
        if g is None:
            continue
        if g.Length() < 60:
            continue
        g = simplify(g, False)
        parts = geom_parts_m(g)
        if parts:
            out["lines"].append((-f.GetField(0), parts))

    # Soundings: shoalest cell in a hierarchy of grids; each spot keeps the first zoom it appears at.
    # Skip spots ashore (OSM land) and implausible drying heights (NONNA carries some land values).
    ashore = land_mask(gt, a.shape)
    depth = np.where(valid & ~ashore & (a < 5), -a, np.inf)
    H, W = a.shape
    chosen = {}
    for size, z in SOUND_CELLS:
        c = int(size / PX)
        hh, ww = (H // c) * c, (W // c) * c
        if hh == 0 or ww == 0:
            continue
        blk = depth[:hh, :ww].reshape(hh // c, c, ww // c, c).transpose(0, 2, 1, 3).reshape(hh // c, ww // c, c * c)
        frac = (blk < 1e9).mean(axis=2)
        k = blk.argmin(axis=2)
        for by, bx in zip(*np.nonzero(frac >= 0.6)):
            iy, ix = divmod(int(k[by, bx]), c)
            key = (by * c + iy, bx * c + ix)
            if size < 200 and depth[key] > 40:      # open deep water needs fewer numbers
                continue
            if key not in chosen:
                chosen[key] = (z, float(depth[key]))
    for (iy, ix), (z, d) in chosen.items():
        if iy >= CHUNK or ix >= CHUNK:       # the overlap row/column belongs to the next chunk
            continue
        mx = gt[0] + (ix + 0.5) * gt[1]
        my = gt[3] + (iy + 0.5) * gt[5]
        out["sound"].append((z, m2w(mx, my), d))
    return out


def land_polygons(box, work, zip_path=None):
    """OSM land polygons clipped to the box, simplified, as lists of world-unit rings."""
    if zip_path is None:
        zip_path = os.path.join(work, "land.zip")
        if not os.path.exists(zip_path):
            print("downloading OSM land polygons", file=sys.stderr)
            get(LAND_URL, zip_path, timeout=1800)
    shp = "/vsizip/" + zip_path + "/land-polygons-split-4326/land_polygons.shp"
    ds = ogr.Open(shp)
    lyr = ds.GetLayer(0)
    w, s, e, n = box
    lyr.SetSpatialFilterRect(w, s, e, n)
    src = osr.SpatialReference(); src.ImportFromEPSG(4326)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst = osr.SpatialReference(); dst.ImportFromEPSG(3857)
    tr = osr.CoordinateTransformation(src, dst)
    clip = ogr.CreateGeometryFromWkt("POLYGON(({w} {s},{e} {s},{e} {n},{w} {n},{w} {s}))".format(w=w, s=s, e=e, n=n))
    out = []
    for f in lyr:
        g = f.GetGeometryRef().Intersection(clip)
        if g is None or g.IsEmpty():
            continue
        g.Transform(tr)
        g = simplify(g)
        parts = geom_parts_m(g)
        if parts:
            out.append(parts)
    return out


# --- OpenStreetMap seamarks -------------------------------------------------------------

COLOUR = {"white": 1, "black": 2, "red": 3, "green": 4, "blue": 5, "yellow": 6, "orange": 11}
SHAPE = {"conical": 1, "can": 2, "spherical": 3, "pillar": 4, "spar": 5, "barrel": 6, "super-buoy": 7}
LAT = {"port": 1, "starboard": 2, "preferred_channel_starboard": 3, "preferred_channel_port": 4}
CARD = {"north": 1, "east": 2, "south": 3, "west": 4}
OBJ = {"buoy_lateral": 1, "buoy_cardinal": 2, "buoy_isolated_danger": 3, "buoy_safe_water": 4,
       "buoy_special_purpose": 5, "beacon_lateral": 6, "beacon_cardinal": 7, "beacon_isolated_danger": 8,
       "beacon_safe_water": 9, "beacon_special_purpose": 10}
LITCHR = {"F": 1, "Fl": 2, "LFl": 3, "Q": 4, "VQ": 5, "UQ": 6, "Iso": 7, "Oc": 8, "IQ": 9, "Mo": 12, "Al": 28}
COL_ABBR = {"white": "W", "red": "R", "green": "G", "yellow": "Y", "blue": "Bu", "orange": "Or"}


def overpass(query):
    last = None
    for url in OVERPASS:
        try:
            data = urllib.parse.urlencode({"data": query}).encode()
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "TillerTime chart builder"})
            with urllib.request.urlopen(req, timeout=900) as r:
                return json.loads(r.read())
        except Exception as e:
            last = e
            print("overpass failed at", url, e, file=sys.stderr)
    raise last


def tag(t, *keys):
    for k in keys:
        if k in t:
            return t[k]
    return None


def light_text(t, pre):
    ch = tag(t, pre + "character") or ""
    grp = tag(t, pre + "group")
    col = tag(t, pre + "colour") or ""
    per = tag(t, pre + "period")
    s = ch + ("(%s)" % grp if grp and grp != "1" else "")
    cols = "".join(COL_ABBR.get(c, "") for c in col.split(";"))
    if cols:
        s += " " + cols
    if per:
        s += " %ss" % per
    return s.strip()


def seamarks(w, box):
    s, wl, n, e = box[1], box[0], box[3], box[2]
    bb = "(%f,%f,%f,%f)" % (s, wl, n, e)
    q = ("[out:json][timeout:800];(node[\"seamark:type\"]%s;way[\"man_made\"~\"^(pier|breakwater|groyne)$\"]%s;"
         "way[\"waterway\"=\"dock\"]%s;);out geom;" % (bb, bb, bb))
    try:
        data = overpass(q)
    except Exception as e:
        print("seamarks skipped:", e, file=sys.stderr)
        return
    for el in data.get("elements", []):
        t = el.get("tags", {})
        if el["type"] == "way":
            pts = [m2w(*lonlat_to_m(p["lon"], p["lat"])) for p in el.get("geometry", [])]
            if len(pts) < 2:
                continue
            closed = pts[0] == pts[-1] and len(pts) > 3
            if closed and t.get("man_made") == "pier" or t.get("waterway") == "dock":
                w.recw(K_STRUCT, 12, 0, None, None, "", [pts])
            else:
                w.recw(K_STRUCTLINE, 12, 0, None, None, "", [pts])
            continue
        st = t.get("seamark:type", "")
        p = m2w(*lonlat_to_m(el["lon"], el["lat"]))
        name = t.get("seamark:name") or t.get("name") or ""
        if st in OBJ:
            pre = "seamark:" + st + ":"
            cols = [COLOUR.get(c, 0) for c in (t.get(pre + "colour") or "").split(";")][:3] + [0, 0, 0]
            cat = LAT.get(t.get(pre + "category", ""), 0) or CARD.get(t.get(pre + "category", ""), 0)
            shp = SHAPE.get(t.get(pre + "shape", ""), 0)
            if st.startswith("beacon") and shp == 0:
                shp = 0
            attr = OBJ[st] | (cat << 4) | (min(shp, 15) << 8) | (cols[0] << 12) | (cols[1] << 16) | (cols[2] << 20)
            m = re.search(r'"?([0-9]+[A-Z]?|[A-Z]{1,2}[0-9]*)"?$', name.strip())
            w.recw(K_NAVAID, 11, attr, None, None, m.group(1) if m else "", [[p]])
        if st in ("light_minor", "light_major", "light_vessel") or "seamark:light:colour" in t or "seamark:light:1:colour" in t:
            pre = "seamark:light:" if "seamark:light:colour" in t else "seamark:light:1:"
            col = COLOUR.get((t.get(pre + "colour") or "white").split(";")[0], 1)
            chv = LITCHR.get(t.get(pre + "character", ""), 0)
            w.recw(K_LIGHT, 10 if st == "light_major" else 12, col | (chv << 4), None, None, light_text(t, pre), [[p]])
        if st in ("rock", "wreck", "obstruction"):
            typ = {"rock": 1, "wreck": 2, "obstruction": 3}[st]
            lvl = t.get("seamark:%s:water_level" % st, "")
            watlev = {"awash": 5, "covers": 4, "submerged": 3, "dry": 2, "always_dry": 2}.get(lvl, 3)
            depth = t.get("seamark:%s:height" % st) if False else t.get("depth") or t.get("seamark:%s:depth" % st)
            try:
                dv = float(depth) if depth else None
            except ValueError:
                dv = None
            w.recw(K_HAZARD, 12, typ | (watlev << 4), dv, None, "", [[p]])


def names(w, box):
    bb = "(%f,%f,%f,%f)" % (box[1], box[0], box[3], box[2])
    q = ("[out:json][timeout:600];(nwr[\"place\"~\"^(island|islet)$\"][\"name\"]%s;"
         "nwr[\"natural\"~\"^(bay|strait|cape)$\"][\"name\"]%s;nwr[\"place\"~\"^(city|town|village|hamlet)$\"][\"name\"]%s;);"
         "out center tags;" % (bb, bb, bb))
    try:
        data = overpass(q)
    except Exception as e:
        print("names skipped:", e, file=sys.stderr)
        return
    for el in data.get("elements", []):
        t = el.get("tags", {})
        c = el.get("center") or ({"lat": el["lat"], "lon": el["lon"]} if "lat" in el else None)
        if not c:
            continue
        p = m2w(*lonlat_to_m(c["lon"], c["lat"]))
        pl, nat = t.get("place"), t.get("natural")
        if nat in ("bay", "strait"):
            sea, z = 1, (9 if nat == "strait" else 12)
        elif pl == "island":
            sea, z = 2, 10
        elif pl == "islet":
            sea, z = 2, 13
        elif pl in ("city", "town"):
            sea, z = 2, 9
        else:
            sea, z = 2, 12
        w.recw(K_NAME, z, sea, None, None, t["name"], [[p]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bbox", default="-125.0,48.2,-122.75,50.0", help="west,south,east,north")
    ap.add_argument("--work", default="bc_work")
    ap.add_argument("--land-zip", default=None)
    ap.add_argument("--output", default="charts_bc.ttenc")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--no-osm", action="store_true")
    ap.add_argument("--depths-only", action="store_true", help="stage 1: write only the NONNA depth records")
    ap.add_argument("--depths-from", default=None, help="stage 2: start from a stage-1 file instead of NONNA")
    a = ap.parse_args()
    box = tuple(map(float, a.bbox.split(",")))
    os.makedirs(a.work, exist_ok=True)
    mx0, my0 = lonlat_to_m(box[0], box[1])
    mx1, my1 = lonlat_to_m(box[2], box[3])
    step = CHUNK * PX
    jobs = []
    y = my1
    while y > my0:
        x = mx0
        while x < mx1:
            jobs.append((x, y, a.work))
            x += step
        y -= step
    print("chunks:", len(jobs), file=sys.stderr)

    w = MWriter()
    if a.depths_from:
        d = open(a.depths_from, "rb").read()
        hl = struct.unpack("<I", d[4:8])[0]
        h = json.loads(d[8:8 + hl])
        w.buf += zlib.decompress(d[8 + hl:], -15)
        w.count = {int(k): v for k, v in h["counts"].items()}
    else:
        load_land_layer(box, a.land_zip or os.path.join(a.work, "land.zip"))
        # Base water tint for the whole region; depth tints and land paint over it.
        corners = [m2w(mx0, my1), m2w(mx1, my1), m2w(mx1, my0), m2w(mx0, my0)]
        w.recw(K_COVER, 0, 0, None, None, "", [corners + [corners[0]]])
        w.recw(K_DEPARE, 0, 0, None, None, "", [corners + [corners[0]]])

        sounds = []
        done = 0
        # The app paints all depth tints, then land, then lines, so record order does not matter.
        with cf.ThreadPoolExecutor(a.threads) as ex:
            for r in ex.map(process_chunk, jobs):
                done += 1
                if done % 20 == 0:
                    print("  %d/%d" % (done, len(jobs)), file=sys.stderr)
                if not r:
                    continue
                for d1, parts in r["areas"]:
                    w.recw(K_DEPARE, 0, 0, d1, None, "", parts)
                for d, parts in r["lines"]:
                    w.recw(K_DEPCNT, 10 if d >= 10 else 11, 0, d, None, "", parts)
                sounds += r["sound"]
        by_zoom = {}
        for z, p, d in sounds:
            by_zoom.setdefault(z, ([], []))
            by_zoom[z][0].append(p); by_zoom[z][1].append(d)
        for z, (pts, ds) in sorted(by_zoom.items()):
            w.recw(K_SOUND, z, 0, None, None, "", [pts], ds)

    if a.depths_only:
        write(w, box, a.output, "CHS NONNA-10 (OGL-Canada)")
        return
    print("land", file=sys.stderr)
    for parts in land_polygons(box, a.work, a.land_zip):
        w.recw(K_LAND, 0, 0, None, None, "", parts)
    if not a.no_osm:
        seamarks(w, box)
        names(w, box)
    write(w, box, a.output, "CHS NONNA-10 (OGL-Canada); OpenStreetMap (ODbL)")


def write(w, box, output, source):
    body = zlib.compress(bytes(w.buf), 9)[2:-4]
    header = json.dumps({"magic": "HENC", "version": 2, "units": "merc30", "side": "ca",
                         "source": source,
                         "generated": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
                         "bbox": box, "records": sum(w.count.values()),
                         "counts": {str(k): v for k, v in sorted(w.count.items())}, "rawBytes": len(w.buf)}).encode()
    with open(output, "wb") as f:
        f.write(b"HENC" + struct.pack("<I", len(header)) + header + body)
    print("records", w.count, "raw %.1f MB, packed %.1f MB" % (len(w.buf) / 1e6, os.path.getsize(output) / 1e6),
          file=sys.stderr)


if __name__ == "__main__":
    main()
