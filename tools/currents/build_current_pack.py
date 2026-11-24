#!/usr/bin/env python3
"""
Build a Tiller Time current pack (.ttcur) from the SalishSeaCast (UBC) forecast.

Source: SalishSeaCast NEMO model, dataset ubcSSfDepthAvgdCurrents1h on
https://salishsea.eos.ubc.ca/erddap/ (Apache License 2.0; cite Soontiens et al. 2016,
Soontiens and Allen 2017). Velocities are already rotated to east/north and averaged
over the upper 5 m, hourly.

What it does
  1. Downloads the model grid (lat/lon of every cell) once and caches it.
  2. Downloads the hourly forecast from `--hours-back` hours ago to the end of the run.
  3. Linearly interpolates the curvilinear model grid onto a regular lat/lon grid
     (~450 m), never across land: triangles with an edge longer than 1.2 km are dropped.
  4. Writes a compact pack: JSON header + one raw-DEFLATE block of
       int32 index map (regular cell -> water cell index, -1 = land)
       then per hour: int16 u[nWater], int16 v[nWater]   (cm/s)

Requires: Python 3.9+, numpy, scipy  (pip install numpy scipy)
Usage:    python3 build_current_pack.py -o salishsea.ttcur
"""
import argparse
import datetime as dt
import json
import os
import struct
import sys
import tempfile
import time
import urllib.request
import zlib

import numpy as np
from scipy.io import netcdf_file
from scipy.spatial import Delaunay

ERDDAP = "https://salishsea.eos.ubc.ca/erddap/griddap"
GRID_DS = "ubcSSnBathymetryV21-08"
CUR_DS = "ubcSSfDepthAvgdCurrents1h"
NY, NX = 898, 398
UA = "TillerTimeCurrentPack/0.1 (+chartplotter; contact via GitHub)"


def fetch(url, path, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=300) as r, open(path, "wb") as f:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b)
            return
        except Exception as e:  # network hiccup: back off and retry
            if k == tries - 1:
                raise
            print(f"  retry after error: {e}", file=sys.stderr)
            time.sleep(5 * (k + 1))


def q(s):
    return s.replace("[", "%5B").replace("]", "%5D")


def load_grid(cache_dir):
    path = os.path.join(cache_dir, "ssc_grid.nc")
    if not os.path.exists(path):
        print("Downloading model grid…")
        url = q(f"{ERDDAP}/{GRID_DS}.nc?latitude[0:1:{NY-1}][0:1:{NX-1}],"
                f"longitude[0:1:{NY-1}][0:1:{NX-1}],bathymetry[0:1:{NY-1}][0:1:{NX-1}]")
        fetch(url, path)
    f = netcdf_file(path, "r", mmap=False)
    lat = np.array(f.variables["latitude"][:], dtype=np.float64)
    lon = np.array(f.variables["longitude"][:], dtype=np.float64)
    f.close()
    return lat, lon


def dataset_end():
    """Last time step currently available (UTC datetime)."""
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "t.csv")
        fetch(q(f"{ERDDAP}/{CUR_DS}.csv?time[(last)]"), p)
        line = open(p).read().strip().splitlines()[-1]
    return dt.datetime.strptime(line.strip(), "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)


def load_currents(t0, t1, workdir, chunk_hours=12):
    """Download VelEast5/VelNorth5 between t0 and t1 inclusive, in chunks."""
    times, us, vs = [], [], []
    cur = t0
    n = 0
    while cur <= t1:
        end = min(cur + dt.timedelta(hours=chunk_hours - 1), t1)
        a, b = cur.strftime("%Y-%m-%dT%H:%M:%SZ"), end.strftime("%Y-%m-%dT%H:%M:%SZ")
        sel = f"[({a}):1:({b})][0:1:{NY-1}][0:1:{NX-1}]"
        url = q(f"{ERDDAP}/{CUR_DS}.nc?VelEast5{sel},VelNorth5{sel}")
        path = os.path.join(workdir, f"cur{n}.nc")
        print(f"Downloading currents {a} → {b}")
        fetch(url, path)
        f = netcdf_file(path, "r", mmap=False)
        times.extend(float(t) for t in f.variables["time"][:])
        us.append(np.array(f.variables["VelEast5"][:], dtype=np.float32))
        vs.append(np.array(f.variables["VelNorth5"][:], dtype=np.float32))
        f.close()
        os.remove(path)
        cur = end + dt.timedelta(hours=1)
        n += 1
    return np.array(times), np.concatenate(us), np.concatenate(vs)


def build(args):
    cache = args.cache or os.path.join(os.path.dirname(os.path.abspath(__file__)), ".cache")
    os.makedirs(cache, exist_ok=True)
    lat, lon = load_grid(cache)

    end = dataset_end()
    now = dt.datetime.now(dt.timezone.utc)
    start = (now - dt.timedelta(hours=args.hours_back)).replace(minute=30, second=0, microsecond=0)
    if start > end:
        start = end - dt.timedelta(hours=24)
    with tempfile.TemporaryDirectory() as work:
        times, U, V = load_currents(start, end, work)
    print(f"{len(times)} hourly frames")

    # Model cells that have data in every frame are water.
    wet = np.isfinite(U).all(axis=0) & np.isfinite(V).all(axis=0)
    src_lat, src_lon = lat[wet], lon[wet]
    print(f"{wet.sum()} wet model cells")

    # Regular target grid
    lat0, lat1 = float(lat.min()), float(lat.max())
    lon0, lon1 = float(lon.min()), float(lon.max())
    dlat, dlon = args.dlat, args.dlon
    nlat = int(np.ceil((lat1 - lat0) / dlat)) + 1
    nlon = int(np.ceil((lon1 - lon0) / dlon)) + 1

    # Local planar coords (km) for triangulation
    c = np.cos(np.radians(0.5 * (lat0 + lat1)))
    def xy(la, lo):
        return np.column_stack(((lo - lon0) * 111.32 * c, (la - lat0) * 110.57))
    P = xy(src_lat, src_lon)
    tri = Delaunay(P)
    S = tri.simplices
    e = lambda i, j: np.hypot(*(P[S[:, i]] - P[S[:, j]]).T)
    longest = np.maximum(np.maximum(e(0, 1), e(1, 2)), e(2, 0))
    good_tri = longest <= args.max_edge_km

    glat = lat0 + dlat * np.arange(nlat)
    glon = lon0 + dlon * np.arange(nlon)
    GLA, GLO = np.meshgrid(glat, glon, indexing="ij")
    T = xy(GLA.ravel(), GLO.ravel())
    simplex = tri.find_simplex(T)
    ok = simplex >= 0
    ok[ok] = good_tri[simplex[ok]]
    # barycentric weights
    idx = np.nonzero(ok)[0]
    s = simplex[idx]
    Tm = tri.transform[s]
    b2 = np.einsum("ijk,ik->ij", Tm[:, :2, :], T[idx] - Tm[:, 2, :])
    w = np.column_stack((b2, 1 - b2.sum(axis=1))).astype(np.float32)
    verts = S[s]
    print(f"regular grid {nlat}×{nlon}, {len(idx)} water cells")

    index_map = np.full(nlat * nlon, -1, dtype=np.int32)
    index_map[idx] = np.arange(len(idx), dtype=np.int32)

    Uw = U[:, wet]
    Vw = V[:, wet]
    nT = len(times)
    out_u = np.empty((nT, len(idx)), dtype=np.int16)
    out_v = np.empty((nT, len(idx)), dtype=np.int16)
    for t in range(nT):
        ut = (Uw[t][verts] * w).sum(axis=1)
        vt = (Vw[t][verts] * w).sum(axis=1)
        out_u[t] = np.clip(np.round(ut * 100), -32767, 32767).astype(np.int16)
        out_v[t] = np.clip(np.round(vt * 100), -32767, 32767).astype(np.int16)

    raw = bytearray(index_map.astype("<i4").tobytes())
    for t in range(nT):
        raw += out_u[t].astype("<i2").tobytes()
        raw += out_v[t].astype("<i2").tobytes()
    comp = zlib.compressobj(9, zlib.DEFLATED, -15)   # raw DEFLATE (Apple COMPRESSION_ZLIB)
    body = comp.compress(bytes(raw)) + comp.flush()

    header = {
        "format": "ttcur",
        "version": 1,
        "source": "SalishSeaCast NEMO (UBC EOAS) " + CUR_DS,
        "variable": "near-surface current, upper 5 m average",
        "license": "Apache-2.0; SalishSeaCast Project Contributors and The University of British Columbia",
        "generated": now.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "lat0": lat0, "lon0": lon0, "dlat": dlat, "dlon": dlon,
        "nlat": nlat, "nlon": nlon, "nWater": int(len(idx)),
        "times": [int(t) for t in times],          # epoch seconds, UTC hour centres
        "scale": 0.01,                              # int16 -> m/s
        "rawBytes": len(raw),
    }
    hj = json.dumps(header, separators=(",", ":")).encode()
    with open(args.output, "wb") as f:
        f.write(b"HCUR")
        f.write(struct.pack("<HHI", 1, 0, len(hj)))
        f.write(hj)
        f.write(body)
    first = dt.datetime.fromtimestamp(times[0], dt.timezone.utc)
    last = dt.datetime.fromtimestamp(times[-1], dt.timezone.utc)
    print(f"Wrote {args.output}: {os.path.getsize(args.output)/1e6:.1f} MB, "
          f"{first:%Y-%m-%d %H:%M} → {last:%Y-%m-%d %H:%M} UTC")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="salishsea.ttcur")
    ap.add_argument("--hours-back", type=int, default=3, help="start this many hours before now")
    ap.add_argument("--dlat", type=float, default=0.004, help="target grid spacing, degrees latitude")
    ap.add_argument("--dlon", type=float, default=0.006, help="target grid spacing, degrees longitude")
    ap.add_argument("--max-edge-km", type=float, default=1.2, help="never interpolate across gaps wider than this")
    ap.add_argument("--cache", help="cache directory for the model grid")
    build(ap.parse_args())
