#!/usr/bin/env python3
"""
Build a Tiller Time tide pack (salishsea.tides.json): high and low tide predictions for
every tide station in the Salish Sea area, US and Canada, for the next 8 days.

Sources (both free, no key):
  NOAA CO-OPS   api.tidesandcurrents.noaa.gov   predictions, interval=hilo, datum MLLW
  CHS IWLS      api-iwls.dfo-mpo.gc.ca          time series wlp-hilo, chart datum
Heights are stored in metres above each station's chart datum (MLLW in the US, Lower Low
Water Large Tide in Canada). The app draws the curve between highs and lows with the
standard cosine shape.

Usage:  python3 build_tide_pack.py -o salishsea.tides.json
"""
import argparse
import datetime as dt
import json
import sys
import time
import urllib.parse
import urllib.request

BOX = dict(south=46.8, north=51.2, west=-126.5, east=-121.2)
UA = "TillerTimeTidePack/0.1 (+chartplotter)"
DAYS = 8


def get_json(url, tries=3):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


def inside(lat, lon):
    return BOX["south"] < lat < BOX["north"] and BOX["west"] < lon < BOX["east"]


def noaa(start):
    out = []
    meta = get_json("https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations.json?type=tidepredictions")
    stations = [s for s in meta["stations"] if inside(s["lat"], s["lng"])]
    print(f"NOAA: {len(stations)} stations")
    for s in stations:
        q = urllib.parse.urlencode(dict(
            product="predictions", application="TillerTime", begin_date=start.strftime("%Y%m%d"),
            range=24 * DAYS, datum="MLLW", station=s["id"], time_zone="gmt", units="metric",
            interval="hilo", format="json"))
        try:
            d = get_json("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?" + q)
        except Exception as e:
            print(f"  skip {s['id']} {s['name']}: {e}", file=sys.stderr)
            continue
        preds = d.get("predictions") or []
        events = []
        for p in preds:
            t = dt.datetime.strptime(p["t"], "%Y-%m-%d %H:%M").replace(tzinfo=dt.timezone.utc)
            events.append([int(t.timestamp()), round(float(p["v"]), 3), "H" if p["type"] == "H" else "L"])
        if len(events) >= 4:
            out.append(dict(id="noaa:" + s["id"], name=s["name"], src="NOAA", lat=s["lat"], lon=s["lng"],
                            datum="MLLW", primary=s.get("type") == "R", events=events))
        time.sleep(0.05)
    return out


def chs(start):
    out = []
    stations = get_json("https://api-iwls.dfo-mpo.gc.ca/api/v1/stations?chs-region-code=PAC")
    stations = [s for s in stations if inside(s["latitude"], s["longitude"])
                and any(ts["code"] == "wlp-hilo" for ts in s.get("timeSeries", []))]
    print(f"CHS: {len(stations)} stations")
    frm = start.strftime("%Y-%m-%dT%H:%M:%SZ")
    to = (start + dt.timedelta(days=DAYS)).strftime("%Y-%m-%dT%H:%M:%SZ")
    for s in stations:
        url = (f"https://api-iwls.dfo-mpo.gc.ca/api/v1/stations/{s['id']}/data?"
               + urllib.parse.urlencode({"time-series-code": "wlp-hilo", "from": frm, "to": to}))
        try:
            d = get_json(url)
        except Exception as e:
            print(f"  skip {s['code']} {s['officialName']}: {e}", file=sys.stderr)
            continue
        pts = []
        for p in d:
            t = dt.datetime.strptime(p["eventDate"][:19], "%Y-%m-%dT%H:%M:%S").replace(tzinfo=dt.timezone.utc)
            pts.append((int(t.timestamp()), float(p["value"])))
        pts.sort()
        events = []
        for i, (t, v) in enumerate(pts):          # label highs and lows by their neighbours
            prev = pts[i - 1][1] if i > 0 else None
            nxt = pts[i + 1][1] if i + 1 < len(pts) else None
            ref = nxt if nxt is not None else prev
            events.append([t, round(v, 3), "H" if ref is None or v > ref else "L"])
        if len(events) >= 4:
            out.append(dict(id="chs:" + s["code"], name=s["officialName"], src="CHS",
                            lat=s["latitude"], lon=s["longitude"], datum="Chart datum (LLWLT)",
                            primary=s.get("type") == "PERMANENT", events=events))
        time.sleep(0.35)                           # IWLS asks for gentle use
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="salishsea.tides.json")
    args = ap.parse_args()
    start = dt.datetime.now(dt.timezone.utc).replace(hour=0, minute=0, second=0, microsecond=0) - dt.timedelta(days=1)
    stations = noaa(start) + chs(start)
    pack = dict(format="tttides", version=1,
                generated=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                units="m", sources=["NOAA CO-OPS (public domain)",
                                    "Canadian Hydrographic Service IWLS (Open Government Licence – Canada)"],
                stations=stations)
    with open(args.output, "w") as f:
        json.dump(pack, f, separators=(",", ":"))
    n = sum(len(s["events"]) for s in stations)
    print(f"Wrote {args.output}: {len(stations)} stations, {n} highs/lows")


if __name__ == "__main__":
    main()
