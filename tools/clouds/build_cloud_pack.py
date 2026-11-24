#!/usr/bin/env python3
"""
Build Tiller Time cloud and outlook packs from free NOAA models (NOMADS grib filter).

1. salishsea.ttcloud  — HRRR 3 km low / mid / high cloud cover, hourly, ~48 h
       "HCLD" u16 1 u16 0 u32 hlen  header JSON  raw DEFLATE: per hour uint8 low[n], mid[n], high[n]  (%)
2. salishsea.ttoutlook — GFS 0.25° 10-day outlook, 3-hourly
       "HOUT" ... per time int16 low[n], mid[n], high[n] (%), u[n], v[n], gust[n] (cm/s), rain[n] (mm/h × 100)
Both on regular lat/lon grids, row 0 = south.

Data: NOAA/NCEP HRRR and GFS (US public domain).
Requires: numpy, scipy, eccodes   (pip install numpy scipy eccodes)
Usage:    python3 build_cloud_pack.py --clouds salishsea.ttcloud --outlook salishsea.ttoutlook
"""
import argparse
import datetime as dt
import json
import os
import struct
import sys
import tempfile
import time
import urllib.parse
import urllib.request
import zlib

import eccodes
import numpy as np
from scipy.spatial import cKDTree

BOX = dict(south=46.8, north=51.2, west=-126.5, east=-121.2)
NOMADS = "https://nomads.ncep.noaa.gov"
UA = "TillerTimeCloudPack/0.1 (+chartplotter)"


def get(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as r:
                d = r.read()
            if not d.startswith(b"GRIB"):
                raise IOError("not GRIB: " + d[:80].decode(errors="replace"))
            return d
        except Exception as e:
            if k == tries - 1:
                raise
            print(f"  retry: {e}", file=sys.stderr)
            time.sleep(4 * (k + 1))


def exists(url):
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status == 200
    except Exception:
        return False


def latest_cycle(model, last_hour, hours=(18, 12, 6, 0)):
    """Newest run whose last forecast file is already on NOMADS."""
    now = dt.datetime.now(dt.timezone.utc)
    for back in range(0, 3):
        day = now - dt.timedelta(days=back)
        for c in hours:
            run = day.replace(hour=c, minute=0, second=0, microsecond=0)
            if run > now:
                continue
            if model == "gfs":
                u = (f"{NOMADS}/pub/data/nccf/com/gfs/prod/gfs.{run:%Y%m%d}/{c:02d}/atmos/"
                     f"gfs.t{c:02d}z.pgrb2.0p25.f{last_hour:03d}.idx")
            else:
                u = (f"{NOMADS}/pub/data/nccf/com/hrrr/prod/hrrr.{run:%Y%m%d}/conus/"
                     f"hrrr.t{c:02d}z.wrfsfcf{last_hour:02d}.grib2.idx")
            if exists(u):
                return run
    raise SystemExit(f"no complete {model} run found")


def messages(data):
    """Yield (shortName, typeOfLevel, stepRange, lats, lons, values) for each GRIB message."""
    with tempfile.NamedTemporaryFile(suffix=".grib2", delete=False) as f:
        f.write(data)
        path = f.name
    try:
        with open(path, "rb") as fh:
            while True:
                g = eccodes.codes_grib_new_from_file(fh)
                if g is None:
                    break
                try:
                    yield (eccodes.codes_get(g, "shortName"), eccodes.codes_get(g, "typeOfLevel"),
                           str(eccodes.codes_get(g, "stepRange")),
                           np.array(eccodes.codes_get_array(g, "latitudes")),
                           np.array(eccodes.codes_get_array(g, "longitudes")),
                           np.array(eccodes.codes_get_values(g)))
                finally:
                    eccodes.codes_release(g)
    finally:
        os.remove(path)


def write_pack(path, magic, header, raw):
    comp = zlib.compressobj(9, zlib.DEFLATED, -15)
    body = comp.compress(bytes(raw)) + comp.flush()
    header["rawBytes"] = len(raw)
    hj = json.dumps(header, separators=(",", ":")).encode()
    with open(path, "wb") as f:
        f.write(magic)
        f.write(struct.pack("<HHI", 1, 0, len(hj)))
        f.write(hj)
        f.write(body)


CLOUD_LEVELS = {"lowCloudLayer": 0, "middleCloudLayer": 1, "highCloudLayer": 2}
CLOUD_VARS = "var_LCDC=on&var_MCDC=on&var_HCDC=on&lev_low_cloud_layer=on&lev_middle_cloud_layer=on&lev_high_cloud_layer=on"


# ---------------------------------------------------------------- HRRR clouds

def build_clouds(path, hours=48, dlat=0.03, dlon=0.045):
    run = latest_cycle("hrrr", hours)
    print(f"HRRR run {run:%Y-%m-%d %HZ}")
    nlat = int(round((BOX["north"] - BOX["south"]) / dlat)) + 1
    nlon = int(round((BOX["east"] - BOX["west"]) / dlon)) + 1
    glat = BOX["south"] + dlat * np.arange(nlat)
    glon = BOX["west"] + dlon * np.arange(nlon)
    GLA, GLO = np.meshgrid(glat, glon, indexing="ij")
    tree = idx = None
    raw = bytearray()
    times = []
    now = dt.datetime.now(dt.timezone.utc)
    for h in range(hours + 1):
        valid = run + dt.timedelta(hours=h)
        if valid < now - dt.timedelta(hours=3):
            continue
        q = (f"dir=%2Fhrrr.{run:%Y%m%d}%2Fconus&file=hrrr.t{run:%H}z.wrfsfcf{h:02d}.grib2&{CLOUD_VARS}"
             f"&subregion=&toplat={BOX['north'] + 0.3}&leftlon={BOX['west'] - 0.5}"
             f"&rightlon={BOX['east'] + 0.5}&bottomlat={BOX['south'] - 0.3}")
        layers = [None, None, None]
        for name, lev, _, la, lo, vals in messages(get(f"{NOMADS}/cgi-bin/filter_hrrr_2d.pl?{q}")):
            if lev in CLOUD_LEVELS:
                if tree is None:
                    lo = np.where(lo > 180, lo - 360, lo)
                    c = np.cos(np.radians(49))
                    tree = cKDTree(np.column_stack((la, lo * c)))
                    _, idx = tree.query(np.column_stack((GLA.ravel(), GLO.ravel() * c)))
                layers[CLOUD_LEVELS[lev]] = vals[idx]
        if any(l is None for l in layers):
            print(f"  f{h:02d} incomplete, skipped", file=sys.stderr)
            continue
        for l in layers:
            raw += np.clip(np.round(l), 0, 100).astype(np.uint8).tobytes()
        times.append(int(valid.timestamp()))
        time.sleep(0.3)
    header = dict(format="ttcloud", version=1, source=f"NOAA HRRR 3 km run {run:%Y-%m-%dT%HZ}",
                  license="US public domain (NOAA/NCEP)",
                  generated=now.strftime("%Y-%m-%dT%H:%M:%SZ"),
                  lat0=float(glat[0]), lon0=float(glon[0]), dlat=dlat, dlon=dlon, nlat=nlat, nlon=nlon,
                  times=times, fields=["low", "mid", "high"], units="%")
    write_pack(path, b"HCLD", header, raw)
    print(f"Wrote {path}: {os.path.getsize(path) / 1e6:.2f} MB, {len(times)} hours, grid {nlat}×{nlon}")


# ---------------------------------------------------------------- GFS outlook

def build_outlook(path, days=10, step=3):
    last = days * 24
    run = latest_cycle("gfs", last)
    print(f"GFS run {run:%Y-%m-%d %HZ}")
    now = dt.datetime.now(dt.timezone.utc)
    lon_w, lon_e = BOX["west"] % 360, BOX["east"] % 360
    lat = lon = None
    frames = []
    times = []
    prev_avg = None          # (start, end, field) of the last 0-6 style average
    for h in range(step, last + 1, step):
        q = (f"dir=%2Fgfs.{run:%Y%m%d}%2F{run:%H}%2Fatmos&file=gfs.t{run:%H}z.pgrb2.0p25.f{h:03d}"
             f"&{CLOUD_VARS}&var_UGRD=on&var_VGRD=on&lev_10_m_above_ground=on&var_GUST=on&var_PRATE=on"
             f"&lev_surface=on&subregion=&toplat={BOX['north']}&leftlon={lon_w}&rightlon={lon_e}&bottomlat={BOX['south']}")
        f = {}
        for name, lev, rng, la, lo, vals in messages(get(f"{NOMADS}/cgi-bin/filter_gfs_0p25.pl?{q}")):
            if lat is None:
                ulat = np.unique(np.round(la, 4)); ulon = np.unique(np.round(lo, 4))
                lat, lon = ulat, ulon
            grid = vals.reshape(len(np.unique(np.round(la, 4))), -1)
            if la[0] > la[-1]:
                grid = grid[::-1]                    # south → north
            if lev in CLOUD_LEVELS:
                f[("cloud", CLOUD_LEVELS[lev])] = grid
            elif name == "10u":
                f["u"] = grid
            elif name == "10v":
                f["v"] = grid
            elif name == "gust":
                f["gust"] = grid
            elif name == "prate" and "-" in rng:
                a, b = (int(x) for x in rng.split("-"))
                f["prate"] = (a, b, grid)
        need = [("cloud", 0), ("cloud", 1), ("cloud", 2), "u", "v", "gust", "prate"]
        if any(k not in f for k in need):
            print(f"  f{h:03d} incomplete, skipped", file=sys.stderr)
            continue
        a, b, avg = f["prate"]
        if b - a == step or prev_avg is None or prev_avg[0] != a:
            rate = avg                               # already a 3-hour average
        else:                                        # 0-6 style: remove the first 3 hours
            pa, pb, pavg = prev_avg
            rate = (avg * (b - a) - pavg * (pb - pa)) / (b - pb)
        prev_avg = (a, b, avg)
        valid = run + dt.timedelta(hours=h)
        if valid < now - dt.timedelta(hours=3):
            continue
        frames.append([f[("cloud", 0)], f[("cloud", 1)], f[("cloud", 2)],
                       f["u"], f["v"], f["gust"], np.maximum(rate, 0) * 3600])   # kg m-2 s-1 → mm/h
        times.append(int(valid.timestamp()))
        time.sleep(0.4)
    raw = bytearray()
    for fr in frames:
        for i, a in enumerate(fr):
            scale = 1 if i < 3 else 100
            raw += np.clip(np.round(a * scale), -32767, 32767).astype("<i2").tobytes()
    lon_deg = np.where(lon > 180, lon - 360, lon)
    header = dict(format="ttoutlook", version=1, source=f"NOAA GFS 0.25° run {run:%Y-%m-%dT%HZ}",
                  license="US public domain (NOAA/NCEP)", generated=now.strftime("%Y-%m-%dT%H:%M:%SZ"),
                  lat0=float(lat[0]), lon0=float(lon_deg[0]), dlat=float(lat[1] - lat[0]),
                  dlon=float(lon[1] - lon[0]), nlat=int(len(lat)), nlon=int(len(lon)), times=times,
                  fields=["low", "mid", "high", "u", "v", "gust", "rain"],
                  scales=[1, 1, 1, 0.01, 0.01, 0.01, 0.01], units=["%", "%", "%", "m/s", "m/s", "m/s", "mm/h"])
    write_pack(path, b"HOUT", header, raw)
    print(f"Wrote {path}: {os.path.getsize(path) / 1e6:.2f} MB, {len(times)} steps to "
          f"{dt.datetime.fromtimestamp(times[-1], dt.timezone.utc):%a %d %b %HZ}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clouds", help="write the HRRR cloud pack here")
    ap.add_argument("--outlook", help="write the GFS 10-day outlook pack here")
    args = ap.parse_args()
    if args.clouds:
        build_clouds(args.clouds)
    if args.outlook:
        build_outlook(args.outlook)
