#!/usr/bin/env python3
"""
Build a Tiller Time wind pack (.ttwind) from Environment and Climate Change Canada's
HRDPS 2.5 km model, through MSC GeoMet (free; Open Government Licence – Canada).

For every forecast hour from `--hours-back` ago to the end of the latest run (~48 h ahead)
it fetches, over the Salish Sea box:
    HRDPS.CONTINENTAL_WSPD        10 m wind speed, m/s
    HRDPS.CONTINENTAL_WD          10 m wind direction, degrees, direction the wind comes FROM
    HRDPS-WEonG_2.5km_WindGust    gust, m/s
and writes:
    "HWND"  u16 version=1  u16 reserved  u32 headerLength  header JSON
    raw DEFLATE block of, per hour:  int16 u[n], int16 v[n], int16 gust[n]   (× 0.01 = m/s)
where u/v are the east/north components of the air motion (toward), n = nlat*nlon,
row 0 = southernmost latitude.

Requires: Python 3.9+, numpy, scipy
Usage:    python3 build_wind_pack.py -o salishsea.ttwind --rain-output salishsea.ttrain

With --rain-output it also fetches HRDPS.CONTINENTAL_RT (precipitation rate) and writes
    "HRAN" ... header JSON ... raw DEFLATE of, per hour: uint8 rate[n]   (× 0.1 = mm/h)
"""
import argparse
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
from scipy.io import netcdf_file

GEOMET = "https://geo.weather.gc.ca/geomet"
SPEED, DIR, GUST = "HRDPS.CONTINENTAL_WSPD", "HRDPS.CONTINENTAL_WD", "HRDPS-WEonG_2.5km_WindGust"
RATE = "HRDPS.CONTINENTAL_RT"          # precipitation rate, mm/s (kg m-2 s-1)
UA = "TillerTimeWindPack/0.1 (+chartplotter)"
BOX = dict(south=46.8, north=51.2, west=-126.5, east=-121.2)   # same area as the current pack


def get(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:
            if k == tries - 1:
                raise
            print(f"  retry after error: {e}", file=sys.stderr)
            time.sleep(3 * (k + 1))


def time_range(layer):
    """(first, last) forecast valid times of the latest run, from the WMS capabilities."""
    xml = get(f"{GEOMET}?service=WMS&version=1.3.0&request=GetCapabilities&layer={layer}").decode()
    m = re.search(r'<Dimension name="time"[^>]*>([^<]+)</Dimension>', xml)
    first, last, _ = m.group(1).split("/")
    p = lambda s: dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)
    ref = re.search(r'<Dimension name="reference_time"[^>]*default="([^"]+)"', xml)
    return p(first), p(last), ref.group(1) if ref else None


def field(layer, t, work):
    url = (f"{GEOMET}?service=WCS&version=2.0.1&request=GetCoverage&coverageId={layer}"
           f"&subset=lat({BOX['south']},{BOX['north']})&subset=long({BOX['west']},{BOX['east']})"
           f"&format=image/netcdf&TIME={t:%Y-%m-%dT%H:%M:%SZ}")
    path = os.path.join(work, "f.nc")
    with open(path, "wb") as f:
        f.write(get(url))
    nc = netcdf_file(path, "r", mmap=False)
    lat = np.array(nc.variables["lat"][:], dtype=np.float64)
    lon = np.array(nc.variables["lon"][:], dtype=np.float64)
    a = np.array(nc.variables["Band1"][:], dtype=np.float32)
    nc.close()
    a[np.abs(a) > 1e30] = np.nan
    if lat[0] > lat[-1]:                 # store south → north
        lat, a = lat[::-1], a[::-1]
    return lat, lon, a


def build(args):
    s0, s1, ref = time_range(SPEED)
    g0, g1, _ = time_range(GUST)
    first, last = max(s0, g0), min(s1, g1)
    now = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0)
    start = max(first, now - dt.timedelta(hours=args.hours_back))
    hours = int((last - start).total_seconds() // 3600) + 1
    print(f"HRDPS run {ref}: {hours} hours {start:%Y-%m-%d %H:%M} → {last:%Y-%m-%d %H:%M} UTC")

    times, U, V, G, P = [], [], [], [], []
    lat = lon = None
    with tempfile.TemporaryDirectory() as work:
        for h in range(hours):
            t = start + dt.timedelta(hours=h)
            la, lo, spd = field(SPEED, t, work)
            _, _, wd = field(DIR, t, work)
            _, _, gust = field(GUST, t, work)
            _, _, rate = field(RATE, t, work) if args.rain_output else (None, None, None)
            if lat is None:
                lat, lon = la, lo
            elif la.shape != lat.shape or lo.shape != lon.shape:
                raise SystemExit("grid changed between hours")
            # direction FROM → velocity components TOWARD
            rad = np.radians(wd)
            U.append(-spd * np.sin(rad))
            V.append(-spd * np.cos(rad))
            G.append(np.where(np.isfinite(gust), np.maximum(gust, spd), spd))
            times.append(int(t.timestamp()))
            if rate is not None:
                P.append(np.nan_to_num(rate, nan=0) * 3600)          # → mm/h
            if h % 6 == 0:
                print(f"  {t:%a %H:%M}Z  max {np.nanmax(spd) * 1.943844:.0f} kn")

    q = lambda a: np.clip(np.nan_to_num(np.round(a * 100), nan=0), -32767, 32767).astype("<i2").tobytes()
    raw = bytearray()
    for u, v, g in zip(U, V, G):
        raw += q(u) + q(v) + q(g)
    comp = zlib.compressobj(9, zlib.DEFLATED, -15)
    body = comp.compress(bytes(raw)) + comp.flush()

    dlat = float(np.mean(np.diff(lat)))
    dlon = float(np.mean(np.diff(lon)))
    header = {
        "format": "ttwind", "version": 1,
        "source": "ECCC HRDPS 2.5 km via MSC GeoMet (" + SPEED + ", " + DIR + ", " + GUST + ")",
        "run": ref,
        "license": "Open Government Licence - Canada; contains information licensed under it",
        "generated": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "lat0": float(lat[0]), "lon0": float(lon[0]), "dlat": dlat, "dlon": dlon,
        "nlat": int(len(lat)), "nlon": int(len(lon)),
        "times": times, "scale": 0.01, "fields": ["u", "v", "gust"],
        "rawBytes": len(raw),
    }
    hj = json.dumps(header, separators=(",", ":")).encode()
    with open(args.output, "wb") as f:
        f.write(b"HWND")
        f.write(struct.pack("<HHI", 1, 0, len(hj)))
        f.write(hj)
        f.write(body)
    print(f"Wrote {args.output}: {os.path.getsize(args.output) / 1e6:.1f} MB, {len(times)} hours, grid {len(lat)}×{len(lon)}")

    if args.rain_output:
        # Rain pack: per hour uint8 rate in 0.1 mm/h steps (0–25.5 mm/h), same grid.
        rraw = bytearray()
        for p in P:
            rraw += np.clip(np.round(p * 10), 0, 255).astype(np.uint8).tobytes()
        comp = zlib.compressobj(9, zlib.DEFLATED, -15)
        rbody = comp.compress(bytes(rraw)) + comp.flush()
        rh = dict(header, format="ttrain", source="ECCC HRDPS 2.5 km via MSC GeoMet (" + RATE + ")",
                  fields=["rate"], scale=0.1, units="mm/h", rawBytes=len(rraw))
        hj = json.dumps(rh, separators=(",", ":")).encode()
        with open(args.rain_output, "wb") as f:
            f.write(b"HRAN")
            f.write(struct.pack("<HHI", 1, 0, len(hj)))
            f.write(hj)
            f.write(rbody)
        peak = max(float(p.max()) for p in P)
        print(f"Wrote {args.rain_output}: {os.path.getsize(args.rain_output) / 1e6:.2f} MB, peak {peak:.1f} mm/h")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="salishsea.ttwind")
    ap.add_argument("--hours-back", type=int, default=3)
    ap.add_argument("--rain-output", help="also write a precipitation-rate pack (.ttrain)")
    build(ap.parse_args())
