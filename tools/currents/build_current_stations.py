#!/usr/bin/env python3
"""Tidal current predictions at stations along the west coast, ten days ahead.

  * US: NOAA CO-OPS current predictions (public domain): max flood, max ebb and slack times.
  * Canada: CHS IWLS (Open Government Licence - Canada): the same events, with the flood and
    ebb directions taken from the station's predicted direction series.

Writes currents_stations.json:
  {"stations": [{"id", "name", "src", "lat", "lon", "flood", "ebb",
                 "ev": [[unix, knots], ...]}]}
knots are signed: + flood, - ebb, 0 slack; flood/ebb are the directions the water flows
toward (degrees true). The app draws the speed between events with the usual cosine curve.
"""
import argparse, concurrent.futures as cf, datetime as dt, json, sys, time, urllib.parse, urllib.request

UA = {"User-Agent": "TillerTime current stations (+chartplotter)"}
BOX = dict(south=32.0, north=61.6, west=-155.0, east=-116.5)
NOAA_LIST = "https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations.json?type=currentpredictions"
NOAA_DATA = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
CHS = "https://api-iwls.dfo-mpo.gc.ca/api/v1"


def get(url, tries=4, timeout=60):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                return json.loads(r.read())
        except Exception as e:
            if k == tries - 1:
                raise
            time.sleep(3 * (k + 1))


def inside(lat, lon):
    return BOX["south"] <= lat <= BOX["north"] and BOX["west"] <= lon <= BOX["east"]


def noaa(start, days):
    st = get(NOAA_LIST)["stations"]
    best = {}
    for s in st:
        if not inside(s["lat"], s["lng"]):
            continue
        best.setdefault(s["id"], (s, None))     # one entry per station; NOAA picks its default bin
    print(f"NOAA: {len(best)} current stations", file=sys.stderr)

    def one(item):
        s, b = item
        q = dict(product="currents_predictions", station=s["id"], begin_date=start.strftime("%Y%m%d"),
                 range=days * 24, interval="MAX_SLACK", units="english", time_zone="gmt", format="json")
        try:
            d = get(NOAA_DATA + "?" + urllib.parse.urlencode(q))
            cp = d.get("current_predictions", {}) if isinstance(d, dict) else {}
            cp = cp.get("cp", []) if isinstance(cp, dict) else []
            if isinstance(cp, dict):
                cp = [cp]
            cp = [p for p in cp if isinstance(p, dict) and "Time" in p]
        except Exception as e:
            print(f"  {s['id']}: {e}", file=sys.stderr)
            return None
        if len(cp) < 4:
            return None
        ev = []
        for p in cp:
            t = dt.datetime.strptime(p["Time"], "%Y-%m-%d %H:%M").replace(tzinfo=dt.timezone.utc)
            v = 0.0 if p["Type"] == "slack" else float(p["Velocity_Major"])
            ev.append([int(t.timestamp()), round(v, 2)])
        f, e = cp[0].get("meanFloodDir"), cp[0].get("meanEbbDir")
        if f is None or e is None:
            return None
        return dict(id=s["id"], name=s["name"].strip(), src="NOAA", lat=round(s["lat"], 5), lon=round(s["lng"], 5),
                    flood=float(f), ebb=float(e), ev=ev)

    out = []
    with cf.ThreadPoolExecutor(6) as ex:
        for r in ex.map(one, best.values()):
            if r:
                out.append(r)
    return out


def chs(start, days):
    stations = get(f"{CHS}/stations")
    cur = [s for s in stations if inside(s["latitude"], s["longitude"]) and
           any(t["code"] == "wcp1-events" for t in s.get("timeSeries", []))]
    print(f"CHS: {len(cur)} current stations", file=sys.stderr)
    out = []
    for s in cur:
        ev, dirs = [], []
        t = start
        end = start + dt.timedelta(days=days)
        try:
            while t < end:                                    # the API serves up to 7 days a call
                t2 = min(t + dt.timedelta(days=6), end)
                a, b = t.strftime("%Y-%m-%dT%H:%M:%SZ"), t2.strftime("%Y-%m-%dT%H:%M:%SZ")
                for p in get(f"{CHS}/stations/{s['id']}/data?time-series-code=wcp1-events&from={a}&to={b}"):
                    tt = dt.datetime.strptime(p["eventDate"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)
                    q = p.get("qualifier", "")
                    v = 0.0 if q == "SLACK" else (abs(p["value"]) if "FLOOD" in q else -abs(p["value"]))
                    ev.append([int(tt.timestamp()), round(v, 2), q])
                t = t2
                time.sleep(0.3)
            # directions: predicted direction near each flood and ebb maximum (first two days)
            a = start.strftime("%Y-%m-%dT%H:%M:%SZ")
            b = (start + dt.timedelta(days=2)).strftime("%Y-%m-%dT%H:%M:%SZ")
            dirs = {int(dt.datetime.strptime(p["eventDate"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)
                        .timestamp()): p["value"]
                    for p in get(f"{CHS}/stations/{s['id']}/data?time-series-code=wcdp1&from={a}&to={b}")}
        except Exception as e:
            print(f"  {s['officialName']}: {e}", file=sys.stderr)
            continue
        if len(ev) < 4 or not dirs:
            continue

        def dir_near(ts):
            k = min(dirs, key=lambda x: abs(x - ts))
            return dirs[k] if abs(k - ts) < 3600 else None

        fl = [dir_near(t) for t, v, q in ev if v > 0]
        eb = [dir_near(t) for t, v, q in ev if v < 0]
        fl = [x for x in fl if x is not None]
        eb = [x for x in eb if x is not None]
        if not fl or not eb:
            continue
        ev.sort()
        out.append(dict(id="CHS-" + s["code"], name=s["officialName"].strip(), src="CHS",
                        lat=round(s["latitude"], 5), lon=round(s["longitude"], 5),
                        flood=float(sorted(fl)[len(fl) // 2]), ebb=float(sorted(eb)[len(eb) // 2]),
                        ev=[[t, v] for t, v, _ in ev]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default="currents_stations.json")
    ap.add_argument("--days", type=int, default=11)
    a = ap.parse_args()
    start = (dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=1)).replace(hour=0, minute=0, second=0, microsecond=0)
    st = noaa(start, a.days) + chs(start, a.days)
    doc = dict(format="ttcurstations", version=1,
               generated=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
               source="NOAA CO-OPS current predictions (public domain); CHS IWLS (Open Government Licence - Canada)",
               stations=st)
    json.dump(doc, open(a.output, "w"), separators=(",", ":"))
    print(f"Wrote {a.output}: {len(st)} stations", file=sys.stderr)


if __name__ == "__main__":
    main()
