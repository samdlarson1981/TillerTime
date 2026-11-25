#!/usr/bin/env python3
"""Build Tiller Time's overview chart pack (HENC v2, band 7): the whole BC and Washington coast
in the same style as the detailed charts, for zoomed-out views and gaps between chart packs.

  * Depth tints and contours: NOAA ETOPO 2022 1 arc-minute (public domain) via ERDDAP.
  * Land: OpenStreetMap land polygons (osmdata.openstreetmap.de), ODbL, simplified to ~40 m.

The app draws it on both sides of the border, underneath the detailed packs, at every zoom.
"""
import argparse, datetime, json, os, struct, sys, urllib.request, zlib

import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
sys.path.insert(0, os.path.dirname(__file__))
import build_bc_pack as bc
from build_enc_pack import K_COVER, K_DEPARE, K_LAND, K_DEPCNT

ETOPO = ("https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.nc?"
         "altitude%5B({s}):1:({n})%5D%5B({w}):1:({e})%5D")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bbox", default="-136.5,45.5,-120.0,56.0")
    ap.add_argument("--work", default="ov_work")
    ap.add_argument("--land-zip", default=None)
    ap.add_argument("--output", default="charts_overview.ttenc")
    a = ap.parse_args()
    box = tuple(map(float, a.bbox.split(",")))
    os.makedirs(a.work, exist_ok=True)
    bc.BAND = 7
    w = bc.MWriter()

    mx0, my0 = bc.lonlat_to_m(box[0], box[1])
    mx1, my1 = bc.lonlat_to_m(box[2], box[3])
    corners = [bc.m2w(mx0, my1), bc.m2w(mx1, my1), bc.m2w(mx1, my0), bc.m2w(mx0, my0)]
    w.recw(K_COVER, 0, 0, None, None, "", [corners + [corners[0]]])
    w.recw(K_DEPARE, 0, 0, None, None, "", [corners + [corners[0]]])

    # ETOPO grid -> Web Mercator raster at ~500 m -> tints and contours.
    nc = os.path.join(a.work, "etopo.nc")
    if not os.path.exists(nc):
        bc.get(ETOPO.format(s=box[1], n=box[3], w=box[0], e=box[2]), nc)
    src = gdal.Open('NETCDF:"%s":altitude' % nc)
    warped = gdal.Warp("", src, format="MEM", dstSRS="EPSG:3857", xRes=500, yRes=500, resampleAlg="cubic",
                       srcSRS="EPSG:4326")
    band = warped.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float64)
    nod = band.GetNoDataValue()
    valid = np.isfinite(arr) if nod is None else (arr != nod) & np.isfinite(arr)
    arr[~valid] = -99999
    band.WriteArray(arr.astype(np.float32))
    band.SetNoDataValue(-99999)
    bc.SIMPLIFY = 120
    s_, lyr = bc.contour(warped, bc.AREA_LEVELS, True)
    for f in lyr:
        g = f.GetGeometryRef()
        emax = f.GetField(1)
        d1 = -1.0 if emax > 1e6 else -emax
        g = bc.simplify(g)
        parts = bc.geom_parts_m(g)
        if parts:
            w.recw(K_DEPARE, 0, 0, d1, None, "", parts)
    s_, lyr = bc.contour(warped, [-200, -100, -50], False)
    for f in lyr:
        g = bc.simplify(f.GetGeometryRef(), False)
        parts = bc.geom_parts_m(g)
        if parts:
            w.recw(K_DEPCNT, 0, 0, -f.GetField(0), None, "", parts)

    bc.SIMPLIFY = 40
    n = 0
    for parts in bc.land_polygons(box, a.work, a.land_zip):
        w.recw(K_LAND, 0, 0, None, None, "", parts)
        n += 1
    print("land polygons", n, file=sys.stderr)

    body = zlib.compress(bytes(w.buf), 9)[2:-4]
    header = json.dumps({"magic": "HENC", "version": 2, "units": "merc30", "side": "all",
                         "source": "NOAA ETOPO 2022 (public domain); OpenStreetMap (ODbL)",
                         "generated": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
                         "bbox": box, "records": sum(w.count.values()),
                         "counts": {str(k): v for k, v in sorted(w.count.items())}, "rawBytes": len(w.buf)}).encode()
    with open(a.output, "wb") as f:
        f.write(b"HENC" + struct.pack("<I", len(header)) + header + body)
    print("records", w.count, "raw %.1f MB, packed %.1f MB" % (len(w.buf) / 1e6, os.path.getsize(a.output) / 1e6),
          file=sys.stderr)


if __name__ == "__main__":
    main()
