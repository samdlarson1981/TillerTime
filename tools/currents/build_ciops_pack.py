#!/usr/bin/env python3
"""
Build a Tiller Time current pack (.ttcur) from Environment and Climate Change Canada's
CIOPS-West ocean model: surface currents (0.5 m) for the BC coast and SE Alaska, hourly
for 48 hours, on a 0.03° × 0.02° grid (~2 km).

Source: MSC Datamart, https://dd.weather.gc.ca/today/model_ciops/west/2km/{HH}/{hhh}/
  ..._MSC_CIOPS-West_SeaWaterVelocity{X,Y}_DBS-0.5m_LatLon0.03x0.02_PT{hhh}H.nc
Contains information licensed under the Open Government Licence – Canada.

Only coastal water is kept (within --coast-km of land), which is what a boat needs and
keeps the pack to a few MB. The app prefers SalishSeaCast (finer) where it has it and the
tidal-current-station mesh in narrow passes this grid can't resolve.

Same format as build_current_pack.py: "HCUR", u16 1, u16 0, u32 header length, JSON header,
then one raw-DEFLATE block of int32 index map (-1 = land) and, per hour, int16 u[nWater],
int16 v[nWater] in cm/s.

Requires: numpy, scipy, netCDF4
Usage:    python3 build_ciops_pack.py -o ciops_west.ttcur
"""
import argparse
import concurrent.futures as cf
import datetime as dt
import json
import os
import re
import struct
import sys
import tempfile
import time
import urllib.request
import zlib

import numpy as np
from scipy.ndimage import distance_transform_edt

BASE = "https://dd.weather.gc.ca/today/model_ciops/west/2km"
NAME = "{run}_MSC_CIOPS-West_SeaWaterVelocity{c}_DBS-0.5m_LatLon0.03x0.02_PT{h:03d}H.nc"


def get(url, tries=4, timeout=120):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            if k == tries - 1:
                raise
            print(f"retry {url}: {e}", file=sys.stderr)
            time.sleep(3 * (k + 1))


def listing(url):
    try:
        return get(url, tries=2, timeout=60).decode("utf-8", "replace")
    except Exception:  # noqa: BLE001
        return ""


def latest_complete_run(hours):
    """Newest run (HH) whose last forecast hour is already published, with its date stamp."""
    runs = sorted(set(re.findall(r'href="(\d\d)/"', listing(BASE + "/"))), reverse=True)
    for hh in runs:
        last = listing(f"{BASE}/{hh}/{hours:03d}/")
        m = re.search(r"(\d{8}T\d\dZ)_MSC_CIOPS-West_SeaWaterVelocityY_DBS-0\.5m", last)
        if m:
            return hh, m.group(1)
    return None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="ciops_west.ttcur")
    ap.add_argument("--hours", type=int, default=48)
    ap.add_argument("--coast-km", type=float, default=15, help="keep water within this distance of land")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    import netCDF4  # noqa: PLC0415

    hh, run = latest_complete_run(args.hours)
    if not hh:
        sys.exit("no complete CIOPS-West run found")
    run_time = dt.datetime.strptime(run, "%Y%m%dT%HZ").replace(tzinfo=dt.timezone.utc)
    print(f"CIOPS-West run {run}", file=sys.stderr)

    tmp = tempfile.mkdtemp(prefix="ciops")
    jobs = [(h, c) for h in range(1, args.hours + 1) for c in "XY"]

    def fetch(job):
        h, c = job
        name = NAME.format(run=run, c=c, h=h)
        path = os.path.join(tmp, name)
        with open(path, "wb") as f:
            f.write(get(f"{BASE}/{hh}/{h:03d}/{name}"))
        return job, path

    paths = {}
    with cf.ThreadPoolExecutor(args.workers) as ex:
        for job, path in ex.map(fetch, jobs):
            paths[job] = path

    def read(path):
        with netCDF4.Dataset(path) as d:
            var = [v for k, v in d.variables.items() if len(v.dimensions) == 3][0]
            return (np.asarray(d["latitude"][:]), np.asarray(d["longitude"][:]),
                    np.ma.filled(var[0].astype(np.float32), np.nan))

    lat, lon, u0 = read(paths[(1, "X")])
    lon = np.where(lon > 180, lon - 360, lon)
    water = np.isfinite(u0)
    near = distance_transform_edt(water) * 2.2 <= args.coast_km     # ~2.2 km per cell
    keep = water & near
    rows = np.where(keep.any(axis=1))[0]
    cols = np.where(keep.any(axis=0))[0]
    r0, r1, c0, c1 = rows.min(), rows.max() + 1, cols.min(), cols.max() + 1
    keep = keep[r0:r1, c0:c1]
    nlat, nlon = keep.shape
    index_map = np.full(keep.shape, -1, dtype=np.int32)
    idx = np.flatnonzero(keep)
    index_map.flat[idx] = np.arange(len(idx), dtype=np.int32)

    raw = bytearray(index_map.astype("<i4").tobytes())
    times = []
    for h in range(1, args.hours + 1):
        u = read(paths[(h, "X")])[2][r0:r1, c0:c1].ravel()[idx]
        v = read(paths[(h, "Y")])[2][r0:r1, c0:c1].ravel()[idx]
        u = np.nan_to_num(u) * 100
        v = np.nan_to_num(v) * 100
        raw += np.clip(np.round(u), -32000, 32000).astype("<i2").tobytes()
        raw += np.clip(np.round(v), -32000, 32000).astype("<i2").tobytes()
        times.append(int((run_time + dt.timedelta(hours=h)).timestamp()))

    comp = zlib.compressobj(9, zlib.DEFLATED, -15)
    body = comp.compress(bytes(raw)) + comp.flush()
    header = {
        "format": "ttcur",
        "version": 1,
        "source": f"CIOPS-West (ECCC) run {run}",
        "variable": "surface current, 0.5 m",
        "license": "Open Government Licence – Canada; Environment and Climate Change Canada",
        "generated": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "lat0": float(lat[r0]), "lon0": float(lon[c0]),
        "dlat": float(lat[1] - lat[0]), "dlon": float(lon[1] - lon[0]),
        "nlat": int(nlat), "nlon": int(nlon), "nWater": int(len(idx)),
        "times": times,
        "scale": 0.01,
        "rawBytes": len(raw),
    }
    hj = json.dumps(header, separators=(",", ":")).encode()
    with open(args.output, "wb") as f:
        f.write(b"HCUR")
        f.write(struct.pack("<HHI", 1, 0, len(hj)))
        f.write(hj)
        f.write(body)
    print(f"Wrote {args.output}: {os.path.getsize(args.output)/1e6:.1f} MB, {len(idx)} cells, "
          f"{len(times)} hours from {run}", file=sys.stderr)


if __name__ == "__main__":
    main()
