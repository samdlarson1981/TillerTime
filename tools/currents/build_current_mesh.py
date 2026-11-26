#!/usr/bin/env python3
"""
Build the tidal-current mesh (.ttmesh): how every patch of coastal water from Puget Sound to
Glacier Bay relates to the tidal current stations near it, so the app can draw a current
field for any time in the ten days the station predictions cover.

For each ~500 m water cell it stores up to three stations, weighted by distance measured
through the water (so a station across an island doesn't count), plus the local channel
direction, so between stations the flow follows the channel rather than a straight-line
average. The app turns this and the stations' predictions (currents_stations.json) into
velocities on the fly; this file only changes when stations or coastlines do.

Land: OpenStreetMap land polygons (© OpenStreetMap contributors, ODbL),
https://osmdata.openstreetmap.de/download/land-polygons-split-4326.zip

Format: "HMSH", u16 1, u16 0, u32 header length, JSON header (grid, station ids), then one
raw-DEFLATE block of
    int32 index[nlat*nlon]          cell -> mesh cell, -1 = land or too far from any station
    u16 s1[n], s2[n], s3[n]         station numbers (65535 = none)
    u8  w1[n], w2[n], w3[n]         weights, summing to 255 × fade (fades out away from stations)
    u8  axis[n]                     channel direction, 0–255 over 0–180°
    u8  strength[n]                 how much to steer along the channel, 0–255
    u8  flags[n]                    bit 0: a fast pass a model grid can't resolve, prefer the mesh

Requires: numpy, scipy, scikit-image, GDAL command-line tools (gdal_rasterize)
Usage:    python3 build_current_mesh.py --land-zip land.zip --stations currents_stations.json -o currents.ttmesh
"""
import argparse
import json
import math
import os
import struct
import subprocess
import sys
import tempfile
import zlib

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter
from skimage.graph import MCP_Geometric

LAT0, LAT1 = 47.0, 59.6
LON0, LON1 = -137.6, -122.0
DLAT, DLON = 0.005, 0.0075


def land_mask(zip_path, nlat, nlon):
    """Rasterize OSM land onto the grid (row 0 = south). True = land."""
    tmp = tempfile.mkdtemp(prefix="mesh")
    tif = os.path.join(tmp, "land.tif")
    src = f"/vsizip/{os.path.abspath(zip_path)}/land-polygons-split-4326/land_polygons.shp"
    clip = os.path.join(tmp, "land.gpkg")
    subprocess.run(["ogr2ogr", "-q", "-f", "GPKG", "-spat", str(LON0 - 0.1), str(LAT0 - 0.1),
                    str(LON1 + 0.1), str(LAT1 + 0.1), clip, src], check=True)
    subprocess.run(["gdal_rasterize", "-q", "-burn", "1", "-ot", "Byte", "-init", "0",
                    "-te", str(LON0 - DLON / 2), str(LAT0 - DLAT / 2), str(LON0 + (nlon - 0.5) * DLON),
                    str(LAT0 + (nlat - 0.5) * DLAT), "-ts", str(nlon), str(nlat),
                    clip, tif], check=True)
    raw = os.path.join(tmp, "land.bin")
    subprocess.run(["gdal_translate", "-q", "-of", "ENVI", tif, raw], check=True)
    a = np.fromfile(raw, dtype=np.uint8).reshape(nlat, nlon)
    return a[::-1] > 0          # GDAL rows run north to south


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--land-zip", required=True)
    ap.add_argument("--stations", required=True, help="currents_stations.json")
    ap.add_argument("-o", "--output", default="currents.ttmesh")
    ap.add_argument("--reach-km", type=float, default=22, help="how far through the water a station counts")
    ap.add_argument("--full-km", type=float, default=7, help="full strength within this of the nearest station")
    args = ap.parse_args()

    nlat = int(round((LAT1 - LAT0) / DLAT)) + 1
    nlon = int(round((LON1 - LON0) / DLON)) + 1
    land = land_mask(args.land_zip, nlat, nlon)
    water = ~land
    print(f"grid {nlat}x{nlon}, water {water.mean():.0%}", file=sys.stderr)

    # Cell size in km (rows are latitude).
    lat_mid = (LAT0 + LAT1) / 2
    dy = DLAT * 111.2
    lats = LAT0 + np.arange(nlat) * DLAT
    dx_row = DLON * 111.2 * np.cos(np.radians(lats))          # varies with latitude
    dx = float(DLON * 111.2 * math.cos(math.radians(lat_mid)))

    st = json.load(open(args.stations))["stations"]
    st = [s for s in st if LAT0 <= s["lat"] <= LAT1 and LON0 <= s["lon"] <= LON1
          and s.get("flood") is not None and s.get("ebb") is not None]
    print(f"{len(st)} stations in the domain", file=sys.stderr)

    costs = np.where(water, 1.0, np.inf)
    K = 3
    best_d = np.full((K, nlat, nlon), np.inf, dtype=np.float32)
    best_s = np.full((K, nlat, nlon), 65535, dtype=np.uint16)
    reach_r = int(args.reach_km / dy) + 2
    for si, s in enumerate(st):
        r = int(round((s["lat"] - LAT0) / DLAT))
        c = int(round((s["lon"] - LON0) / DLON))
        if not (0 <= r < nlat and 0 <= c < nlon):
            continue
        # Snap to the nearest water cell (stations sit right at the shore sometimes).
        if not water[r, c]:
            win = water[max(r - 4, 0):r + 5, max(c - 4, 0):c + 5]
            ys, xs = np.nonzero(win)
            if len(ys) == 0:
                continue
            k = np.argmin((ys + max(r - 4, 0) - r) ** 2 + (xs + max(c - 4, 0) - c) ** 2)
            r, c = ys[k] + max(r - 4, 0), xs[k] + max(c - 4, 0)
        rr0, rr1 = max(r - reach_r, 0), min(r + reach_r + 1, nlat)
        cr = int(args.reach_km / dx_row[r]) + 2
        cc0, cc1 = max(c - cr, 0), min(c + cr + 1, nlon)
        sub = costs[rr0:rr1, cc0:cc1]
        m = MCP_Geometric(sub, sampling=(dy, float(dx_row[r])))
        cum, _ = m.find_costs([(r - rr0, c - cc0)])
        d = cum.astype(np.float32)
        d[~np.isfinite(d) | (d > args.reach_km)] = np.inf
        # Insert into the three nearest (an insertion sort over the window).
        bd = best_d[:, rr0:rr1, cc0:cc1]
        bs = best_s[:, rr0:rr1, cc0:cc1]
        s_arr = np.full(d.shape, si, dtype=np.uint16)
        for k in range(K):
            closer = d < bd[k]
            if not closer.any():
                break
            nd, ns = np.where(closer, bd[k], d), np.where(closer, bs[k], s_arr)
            bd[k][...] = np.where(closer, d, bd[k])
            bs[k][...] = np.where(closer, s_arr, bs[k])
            d, s_arr = nd, ns
        if si % 100 == 0:
            print(f"  station {si}/{len(st)}", file=sys.stderr)

    covered = water & np.isfinite(best_d[0])
    d1 = best_d[0]
    fade = np.clip((args.reach_km - d1) / (args.reach_km - args.full_km), 0, 1)
    wts = np.where(np.isfinite(best_d), 1.0 / np.maximum(best_d, 0.3) ** 2, 0.0)
    wsum = wts.sum(axis=0)
    wsum[wsum == 0] = 1
    wq = np.round(wts / wsum * 255 * fade).astype(np.uint8)
    covered &= wq.sum(axis=0) > 3

    # Channel direction: across the gradient of distance-to-land (in km), from the smoothed
    # structure tensor; strength from how coherent it is and how narrow the water is.
    dist_km = distance_transform_edt(water, sampling=(dy, dx))
    gy, gx = np.gradient(dist_km, dy, dx)
    sig = 1.5
    jxx = gaussian_filter(gx * gx, sig)
    jyy = gaussian_filter(gy * gy, sig)
    jxy = gaussian_filter(gx * gy, sig)
    grad_dir = 0.5 * np.arctan2(2 * jxy, jxx - jyy)              # dominant gradient angle, from east
    lam = np.sqrt((jxx - jyy) ** 2 + 4 * jxy ** 2)
    coh = np.where(jxx + jyy > 1e-9, lam / (jxx + jyy + 1e-9), 0)
    axis_east = grad_dir + np.pi / 2                             # channel runs across the gradient
    axis_true = np.degrees(np.pi / 2 - axis_east) % 180          # compass bearing, 0–180
    width_km = gaussian_filter(dist_km, 2) * 2
    narrow = np.clip((5.0 - width_km) / 3.5, 0, 1)
    strength = np.clip(coh * narrow, 0, 1)

    # Prefer the mesh over the models in fast, narrow passes.
    peaks = np.array([max(abs(e[1]) for e in s["ev"]) if s["ev"] else 0 for s in st], dtype=np.float32)
    s1 = best_s[0]
    fast = np.zeros_like(covered)
    valid = s1 < 65535
    fast[valid] = (peaks[s1[valid]] >= 3) & (d1[valid] <= 2.5) & (width_km[valid] <= 3)

    index_map = np.full((nlat, nlon), -1, dtype=np.int32)
    idx = np.flatnonzero(covered)
    index_map.flat[idx] = np.arange(len(idx), dtype=np.int32)
    raw = bytearray(index_map.astype("<i4").tobytes())
    for k in range(K):
        raw += best_s[k].ravel()[idx].astype("<u2").tobytes()
    for k in range(K):
        raw += wq[k].ravel()[idx].tobytes()
    raw += np.round(axis_true.ravel()[idx] / 180 * 255).clip(0, 255).astype(np.uint8).tobytes()
    raw += np.round(strength.ravel()[idx] * 255).astype(np.uint8).tobytes()
    raw += fast.ravel()[idx].astype(np.uint8).tobytes()

    comp = zlib.compressobj(9, zlib.DEFLATED, -15)
    body = comp.compress(bytes(raw)) + comp.flush()
    header = {
        "format": "ttmesh", "version": 1,
        "source": "Tidal current stations (NOAA CO-OPS, CHS) on an OpenStreetMap water mask",
        "license": "Land: © OpenStreetMap contributors (ODbL)",
        "lat0": LAT0, "lon0": LON0, "dlat": DLAT, "dlon": DLON, "nlat": nlat, "nlon": nlon,
        "nCells": int(len(idx)), "stations": [s["id"] for s in st], "rawBytes": len(raw),
    }
    hj = json.dumps(header, separators=(",", ":")).encode()
    with open(args.output, "wb") as f:
        f.write(b"HMSH")
        f.write(struct.pack("<HHI", 1, 0, len(hj)))
        f.write(hj)
        f.write(body)
    print(f"Wrote {args.output}: {os.path.getsize(args.output)/1e6:.1f} MB, {len(idx)} cells, "
          f"{int(fast.sum())} fast-pass cells", file=sys.stderr)


if __name__ == "__main__":
    main()
