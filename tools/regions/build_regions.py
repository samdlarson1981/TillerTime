#!/usr/bin/env python3
"""Build Tiller Time's wind, rain, cloud and 10-day outlook packs for the whole west coast,
San Diego to Kodiak, as a set of overlapping regions. The app loads the region under the map.

Per region it writes  <id>.ttwind  <id>.ttrain  <id>.ttcloud  <id>.ttoutlook  and a shared
regions.json listing them. The Salish Sea region keeps its old file names (salishsea.*).

Models (all free):
  * ECCC HRDPS 2.5 km (Open Government Licence - Canada) via MSC GeoMet: wind, gusts and rain
    from northern California to Southeast Alaska.
  * NOAA HRRR 3 km CONUS and Alaska (public domain) from the NOAA Open Data bucket on AWS: wind,
    gusts and rain in California and the Gulf of Alaska; low/mid/high cloud everywhere it reaches.
  * NOAA GFS 0.25 deg (public domain) via NOMADS: the 10-day outlook, and clouds on the central BC
    coast where no HRRR domain reaches.

Requires: numpy, scipy, eccodes
"""
import argparse, concurrent.futures as cf, datetime as dt, json, os, re, struct, sys, tempfile, time, \
    urllib.request, zlib

import eccodes
import numpy as np
from scipy.io import netcdf_file
from scipy.spatial import cKDTree

UA = "TillerTimeRegions/0.1 (+chartplotter)"
GEOMET = "https://geo.weather.gc.ca/geomet"
NOMADS = "https://nomads.ncep.noaa.gov"
HRRR_AWS = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"

# id, name, south, north, west, east, wind model ("hrdps" | "conus" | "ak")
REGIONS = [
    ("socal", "Southern California", 32.3, 35.3, -121.3, -116.9, "conus"),
    ("centralca", "Central California", 34.6, 38.7, -124.3, -120.4, "conus"),
    ("northca", "Northern California", 38.3, 42.4, -125.2, -122.2, "conus"),
    ("oregon", "Oregon", 41.8, 46.9, -125.5, -122.5, "hrdps"),
    ("salishsea", "Salish Sea and Washington", 46.8, 51.2, -126.5, -121.2, "hrdps"),
    ("westvi", "West Vancouver Island and Central Coast", 48.3, 52.6, -130.0, -124.0, "hrdps"),
    ("northbc", "North Coast and Haida Gwaii", 51.8, 55.6, -133.6, -127.5, "hrdps"),
    ("seak", "Southeast Alaska", 54.5, 59.8, -139.5, -130.3, "hrdps"),
    ("gulfak", "Gulf of Alaska", 56.8, 61.4, -154.8, -138.0, "ak"),
]
# Grid spacing (degrees latitude) for HRRR-built wind and for clouds; the big Gulf of Alaska
# region uses ~5 km so its download stays a few MB.
GRID_DLAT = {"gulfak": 0.05}
COAST = dict(south=32.0, north=61.6, west=-155.0, east=-116.5)

HRRR_VARS = {  # idx pattern -> field name
    ":UGRD:10 m above ground:": "u", ":VGRD:10 m above ground:": "v", ":GUST:surface:": "gust",
    ":PRATE:surface:": "prate", ":LCDC:low cloud layer:": "low", ":MCDC:middle cloud layer:": "mid",
    ":HCDC:high cloud layer:": "high",
}


def get(url, headers=None, tries=4, timeout=180):
    h = {"User-Agent": UA}
    h.update(headers or {})
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=timeout) as r:
                return r.read()
        except Exception as e:
            if k == tries - 1:
                raise
            print(f"  retry {url[:90]}: {e}", file=sys.stderr)
            time.sleep(4 * (k + 1))


def exists(url):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA}), timeout=30) as r:
            return r.status == 200
    except Exception:
        return False


def write_pack(path, magic, header, raw):
    comp = zlib.compressobj(9, zlib.DEFLATED, -15)
    body = comp.compress(bytes(raw)) + comp.flush()
    header["rawBytes"] = len(raw)
    hj = json.dumps(header, separators=(",", ":")).encode()
    with open(path, "wb") as f:
        f.write(magic + struct.pack("<HHI", 1, 0, len(hj)) + hj + body)


def grib_fields(data):
    """[(lat, lon, values, keys)] for each message in a GRIB2 byte string."""
    out = []
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
                    keys = {k: eccodes.codes_get(g, k) for k in ("shortName", "typeOfLevel", "gridType")}
                    keys["stepRange"] = str(eccodes.codes_get(g, "stepRange"))
                    for k in ("LoVInDegrees", "Latin1InDegrees", "orientationOfTheGridInDegrees",
                              "resolutionAndComponentFlags"):
                        try:
                            keys[k] = eccodes.codes_get(g, k)
                        except Exception:
                            pass
                    out.append((np.array(eccodes.codes_get_array(g, "latitudes")),
                                np.array(eccodes.codes_get_array(g, "longitudes")),
                                np.array(eccodes.codes_get_values(g)), keys))
                finally:
                    eccodes.codes_release(g)
    finally:
        os.remove(path)
    return out


# ------------------------------------------------------------------------- target grids

class Grid:
    def __init__(self, s, n, w, e, dlat, dlon):
        self.lat = s + dlat * np.arange(int(round((n - s) / dlat)) + 1)
        self.lon = w + dlon * np.arange(int(round((e - w) / dlon)) + 1)
        self.dlat, self.dlon = dlat, dlon
        self.LA, self.LO = np.meshgrid(self.lat, self.lon, indexing="ij")

    def meta(self):
        return dict(lat0=float(self.lat[0]), lon0=float(self.lon[0]), dlat=self.dlat, dlon=self.dlon,
                    nlat=len(self.lat), nlon=len(self.lon))


class Mapper:
    """Nearest source point for each target point, with a 'too far' mask (outside the domain)."""

    def __init__(self, slat, slon, grid, max_km=6.0):
        slon = np.where(slon > 180, slon - 360, slon)
        c = np.cos(np.radians(np.mean(grid.lat)))
        tree = cKDTree(np.column_stack((slat, slon * c)))
        d, self.idx = tree.query(np.column_stack((grid.LA.ravel(), grid.LO.ravel() * c)))
        self.ok = (d * 111.0) <= max_km
        self.shape = grid.LA.shape

    def apply(self, vals):
        out = np.full(self.idx.shape, np.nan)
        out[self.ok] = vals[self.idx[self.ok]]
        return out.reshape(self.shape)


# ------------------------------------------------------------------------------- GFS

def latest_gfs(last_hour):
    now = dt.datetime.now(dt.timezone.utc)
    for back in range(3):
        for c in (18, 12, 6, 0):
            run = (now - dt.timedelta(days=back)).replace(hour=c, minute=0, second=0, microsecond=0)
            if run <= now and exists(f"{NOMADS}/pub/data/nccf/com/gfs/prod/gfs.{run:%Y%m%d}/{c:02d}/atmos/"
                                     f"gfs.t{c:02d}z.pgrb2.0p25.f{last_hour:03d}.idx"):
                return run
    raise SystemExit("no complete GFS run")


def gfs_outlook(days=10, step=3):
    days = int(os.environ.get("OUTLOOK_DAYS", days))
    """GFS frames over the whole coast: (lat, lon, times, frames[name] -> [time][lat][lon])."""
    last = days * 24
    run = latest_gfs(last)
    print(f"GFS run {run:%Y-%m-%d %HZ}")
    now = dt.datetime.now(dt.timezone.utc)
    q0 = ("var_LCDC=on&var_MCDC=on&var_HCDC=on&lev_low_cloud_layer=on&lev_middle_cloud_layer=on"
          "&lev_high_cloud_layer=on&var_UGRD=on&var_VGRD=on&lev_10_m_above_ground=on&var_GUST=on&var_PRATE=on"
          f"&lev_surface=on&subregion=&toplat={COAST['north']}&leftlon={COAST['west'] % 360}"
          f"&rightlon={COAST['east'] % 360}&bottomlat={COAST['south']}")
    lat = lon = None
    times, frames, prev = [], {k: [] for k in ("low", "mid", "high", "u", "v", "gust", "rain")}, None
    for h in range(step, last + 1, step):
        q = f"dir=%2Fgfs.{run:%Y%m%d}%2F{run:%H}%2Fatmos&file=gfs.t{run:%H}z.pgrb2.0p25.f{h:03d}&" + q0
        f = {}
        for la, lo, vals, k in grib_fields(get(f"{NOMADS}/cgi-bin/filter_gfs_0p25.pl?{q}")):
            ulat = np.unique(np.round(la, 4)); ulon = np.unique(np.round(lo, 4))
            grid = vals.reshape(len(ulat), len(ulon))
            if la[0] > la[-1]:
                grid = grid[::-1]
            lat, lon = ulat, np.where(ulon > 180, ulon - 360, ulon)
            name, lev, rng = k["shortName"], k["typeOfLevel"], k["stepRange"]
            if lev == "lowCloudLayer": f["low"] = grid
            elif lev == "middleCloudLayer": f["mid"] = grid
            elif lev == "highCloudLayer": f["high"] = grid
            elif name == "10u": f["u"] = grid
            elif name == "10v": f["v"] = grid
            elif name == "gust": f["gust"] = grid
            elif name == "prate" and "-" in rng:
                a, b = (int(x) for x in rng.split("-"))
                f["prate"] = (a, b, grid)
        if any(k not in f for k in ("low", "mid", "high", "u", "v", "gust", "prate")):
            print(f"  GFS f{h:03d} incomplete, skipped", file=sys.stderr)
            continue
        a, b, avg = f["prate"]
        if b - a == step or prev is None or prev[0] != a:
            rate = avg
        else:
            pa, pb, pavg = prev
            rate = (avg * (b - a) - pavg * (pb - pa)) / (b - pb)
        prev = (a, b, avg)
        valid = run + dt.timedelta(hours=h)
        if valid < now - dt.timedelta(hours=3):
            continue
        for k in ("low", "mid", "high", "u", "v", "gust"):
            frames[k].append(f[k])
        frames["rain"].append(np.maximum(rate, 0) * 3600)
        times.append(int(valid.timestamp()))
        time.sleep(0.3)
    return run, lat, lon, times, {k: np.array(v) for k, v in frames.items()}


def geps_pop(gfs, work):
    """Chance of rain (%) on the GFS outlook grid and times: ECCC's global ensemble (GEPS,
    21 members), probability of at least 1 mm in the 6 hours around each time. Whole coast."""
    run, lat, lon, times, fr = gfs
    layer = "GEPS.DIAG.6_PRMM.ERGE1"
    xml = get(f"{GEOMET}?service=WMS&version=1.3.0&request=GetCapabilities&layer={layer}").decode()
    m = re.search(r'<Dimension name="time"[^>]*>([^<]+)</Dimension>', xml)
    ends = [dt.datetime.strptime(t, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)
            for t in m.group(1).split(",")] if m and "," in m.group(1) else []
    if not ends:
        print("  GEPS: no time list, chance of rain skipped", file=sys.stderr)
        return None
    region = ("coast", "", COAST["south"], COAST["north"], COAST["west"], COAST["east"], "")
    out = np.zeros((len(times), len(lat), len(lon)), np.float32)
    cache = {}
    for k, t in enumerate(times):
        tt = dt.datetime.fromtimestamp(t, dt.timezone.utc) + dt.timedelta(hours=3)   # window centred on t
        e = next((x for x in ends if x >= tt), None)
        if e is None:
            out[k] = np.nan
            continue
        if e not in cache:
            try:
                glat, glon, a = hrdps_field(layer, region, e, work)
                cache[e] = (glat, glon, np.nan_to_num(a, nan=0))
            except Exception as ex:
                print(f"  GEPS {e:%d %H}Z: {ex}", file=sys.stderr)
                cache[e] = None
        if cache[e] is None:
            out[k] = np.nan
            continue
        glat, glon, a = cache[e]
        ii = np.clip(np.searchsorted(glat, lat), 1, len(glat) - 1)
        jj = np.clip(np.searchsorted(glon, lon), 1, len(glon) - 1)
        ii = np.where(np.abs(glat[ii - 1] - lat) < np.abs(glat[ii] - lat), ii - 1, ii)
        jj = np.where(np.abs(glon[jj - 1] - lon) < np.abs(glon[jj] - lon), jj - 1, jj)
        out[k] = a[np.ix_(ii, jj)]
        time.sleep(0.1)
    return np.clip(out, 0, 100)


def write_outlook(path, region, gfs, pop=None):
    run, lat, lon, times, fr = gfs
    _, _, s, n, w, e, _ = region
    i0, i1 = np.searchsorted(lat, s - 0.25), np.searchsorted(lat, n + 0.25)
    j0, j1 = np.searchsorted(lon, w - 0.25), np.searchsorted(lon, e + 0.25)
    raw = bytearray()
    for t in range(len(times)):
        for i, k in enumerate(("low", "mid", "high", "u", "v", "gust", "rain")):
            a = fr[k][t, i0:i1, j0:j1]
            raw += np.clip(np.round(a * (1 if i < 3 else 100)), -32767, 32767).astype("<i2").tobytes()
        if pop is not None:      # 8th field: chance of rain %, -1 where unknown
            raw += np.nan_to_num(np.round(pop[t, i0:i1, j0:j1]), nan=-1).astype("<i2").tobytes()
    header = dict(format="ttoutlook", version=1, source=f"NOAA GFS 0.25° run {run:%Y-%m-%dT%HZ}",
                  license="US public domain (NOAA/NCEP)",
                  generated=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                  lat0=float(lat[i0]), lon0=float(lon[j0]), dlat=float(lat[1] - lat[0]), dlon=float(lon[1] - lon[0]),
                  nlat=int(i1 - i0), nlon=int(j1 - j0), times=times,
                  fields=["low", "mid", "high", "u", "v", "gust", "rain"] + (["pop"] if pop is not None else []),
                  scales=[1, 1, 1, 0.01, 0.01, 0.01, 0.01, 1][:8 if pop is not None else 7],
                  units=["%", "%", "%", "m/s", "m/s", "m/s", "mm/h", "%"][:8 if pop is not None else 7],
                  popSource="ECCC GEPS ensemble: chance of >= 1 mm in 6 h (Open Government Licence - Canada)")
    write_pack(path, b"HOUT", header, raw)


def gfs_clouds_at(gfs, grid, t):
    """GFS low/mid/high at unix time t on a target grid (time-linear, nearest in space)."""
    _, lat, lon, times, fr = gfs
    ts = np.array(times)
    k = int(np.clip(np.searchsorted(ts, t), 1, len(ts) - 1))
    w = float(np.clip((t - ts[k - 1]) / max(ts[k] - ts[k - 1], 1), 0, 1))
    ii = np.clip(np.round((grid.LA - lat[0]) / (lat[1] - lat[0])).astype(int), 0, len(lat) - 1)
    jj = np.clip(np.round((grid.LO - lon[0]) / (lon[1] - lon[0])).astype(int), 0, len(lon) - 1)
    return [fr[c][k - 1][ii, jj] * (1 - w) + fr[c][k][ii, jj] * w for c in ("low", "mid", "high")]


# ------------------------------------------------------------------------------ HRRR

def latest_hrrr(domain, last=48):
    now = dt.datetime.now(dt.timezone.utc)
    sfx = ".ak.grib2" if domain == "alaska" else ".grib2"
    for back in range(3):
        for c in (18, 12, 6, 0):
            run = (now - dt.timedelta(days=back)).replace(hour=c, minute=0, second=0, microsecond=0)
            if run <= now and exists(f"{HRRR_AWS}/hrrr.{run:%Y%m%d}/{domain}/hrrr.t{c:02d}z.wrfsfcf{last:02d}{sfx}.idx"):
                return run
    raise SystemExit(f"no complete HRRR {domain} run")


def hrrr_hour(domain, run, h):
    """{name: values on the native grid}, plus (lat, lon) and earth-relative winds."""
    sfx = ".ak.grib2" if domain == "alaska" else ".grib2"
    url = f"{HRRR_AWS}/hrrr.{run:%Y%m%d}/{domain}/hrrr.t{run:%H}z.wrfsfcf{h:02d}{sfx}"
    idx = get(url + ".idx").decode().splitlines()
    out, lat, lon, keys = {}, None, None, None
    for i, line in enumerate(idx):
        for pat, name in HRRR_VARS.items():
            if pat in line and name not in out:
                a = int(line.split(":")[1])
                b = int(idx[i + 1].split(":")[1]) - 1 if i + 1 < len(idx) else ""
                msgs = grib_fields(get(url, {"Range": f"bytes={a}-{b}"}))
                if msgs:
                    la, lo, vals, k = msgs[0]
                    out[name] = vals
                    if lat is None:
                        lat, lon, keys = la, np.where(lo > 180, lo - 360, lo), k
    if not all(n in out for n in HRRR_VARS.values()):
        return None
    # HRRR winds are relative to the projection grid; turn them to true north/east.
    if keys.get("resolutionAndComponentFlags", 0) & 8:
        if keys["gridType"] == "lambert":
            n, lov = np.sin(np.radians(keys["Latin1InDegrees"])), keys["LoVInDegrees"]
        else:                                     # polar stereographic
            n, lov = 1.0, keys["orientationOfTheGridInDegrees"]
        lov = lov - 360 if lov > 180 else lov
        ang = np.radians(n * (lon - lov))
        u, v = out["u"], out["v"]
        out["u"], out["v"] = np.cos(ang) * u + np.sin(ang) * v, -np.sin(ang) * u + np.cos(ang) * v
    return lat, lon, out


# ----------------------------------------------------------------------------- HRDPS

def hrdps_range(layer):
    xml = get(f"{GEOMET}?service=WMS&version=1.3.0&request=GetCapabilities&layer={layer}").decode()
    m = re.search(r'<Dimension name="time"[^>]*>([^<]+)</Dimension>', xml)
    first, last, _ = m.group(1).split("/")
    p = lambda s: dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)
    ref = re.search(r'<Dimension name="reference_time"[^>]*default="([^"]+)"', xml)
    return p(first), p(last), ref.group(1) if ref else None


def hrdps_field(layer, region, t, work):
    _, _, s, n, w, e, _ = region
    url = (f"{GEOMET}?service=WCS&version=2.0.1&request=GetCoverage&coverageId={layer}"
           f"&subset=lat({s},{n})&subset=long({w},{e})&format=image/netcdf&TIME={t:%Y-%m-%dT%H:%M:%SZ}")
    path = os.path.join(work, f"{region[0]}_{layer}.nc")
    open(path, "wb").write(get(url))
    nc = netcdf_file(path, "r", mmap=False)
    lat = np.array(nc.variables["lat"][:], dtype=np.float64)
    lon = np.array(nc.variables["lon"][:], dtype=np.float64)
    a = np.array(nc.variables["Band1"][:], dtype=np.float32)
    nc.close()
    a[np.abs(a) > 1e30] = np.nan
    if lat[0] > lat[-1]:
        lat, a = lat[::-1], a[::-1]
    return lat, lon, a


def build_hrdps(region, start, hours, ref, work):
    times, U, V, G, P = [], [], [], [], []
    lat = lon = None
    for h in range(hours):
        t = start + dt.timedelta(hours=h)
        try:
            la, lo, spd = hrdps_field("HRDPS.CONTINENTAL_WSPD", region, t, work)
            _, _, wd = hrdps_field("HRDPS.CONTINENTAL_WD", region, t, work)
            _, _, gust = hrdps_field("HRDPS-WEonG_2.5km_WindGust", region, t, work)
            _, _, rate = hrdps_field("HRDPS.CONTINENTAL_RT", region, t, work)
        except Exception as e:
            print(f"  {region[0]} {t:%d %H}Z skipped: {e}", file=sys.stderr)
            continue
        if lat is None:
            lat, lon = la, lo
        elif la.shape != lat.shape:
            continue
        rad = np.radians(wd)
        U.append(-spd * np.sin(rad)); V.append(-spd * np.cos(rad))
        G.append(np.where(np.isfinite(gust), np.maximum(gust, spd), spd))
        P.append(np.nan_to_num(rate, nan=0) * 3600)
        times.append(int(t.timestamp()))
    meta = dict(lat0=float(lat[0]), lon0=float(lon[0]), dlat=float(np.mean(np.diff(lat))),
                dlon=float(np.mean(np.diff(lon))), nlat=len(lat), nlon=len(lon))
    src = f"ECCC HRDPS 2.5 km via MSC GeoMet, run {ref}"
    return meta, times, U, V, G, P, src, "Open Government Licence - Canada; contains information licensed under it"


def write_wind_rain(region, out_dir, meta, times, U, V, G, P, src, lic):
    rid = region[0]
    q = lambda a: np.clip(np.nan_to_num(np.round(a * 100), nan=0), -32767, 32767).astype("<i2").tobytes()
    raw = bytearray()
    for u, v, g in zip(U, V, G):
        raw += q(u) + q(v) + q(g)
    gen = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    header = dict(format="ttwind", version=1, source=src, license=lic, generated=gen, region=rid,
                  times=times, scale=0.01, fields=["u", "v", "gust"], **meta)
    write_pack(os.path.join(out_dir, f"{rid}.ttwind"), b"HWND", header, raw)
    rraw = bytearray()
    for p in P:
        rraw += np.clip(np.round(np.nan_to_num(p) * 10), 0, 255).astype(np.uint8).tobytes()
    rh = dict(header, format="ttrain", fields=["rate"], scale=0.1, units="mm/h")
    write_pack(os.path.join(out_dir, f"{rid}.ttrain"), b"HRAN", rh, rraw)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="packs")
    ap.add_argument("--only", default=None, help="comma-separated region ids")
    ap.add_argument("--hours", type=int, default=48)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    regions = [r for r in REGIONS if not a.only or r[0] in a.only.split(",")]
    now = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0)

    gfs = gfs_outlook()
    with tempfile.TemporaryDirectory() as work:
        try:
            pop = geps_pop(gfs, work)
        except Exception as e:
            print("  chance of rain skipped:", e, file=sys.stderr)
            pop = None
    for r in regions:
        write_outlook(os.path.join(a.out, f"{r[0]}.ttoutlook"), r, gfs, pop)
    print("outlooks written")

    # Cloud grids (all regions) and wind grids (HRRR regions) on regular lat/lon.
    def grid(r, d):
        d = GRID_DLAT.get(r[0], d)
        return Grid(r[2], r[3], r[4], r[5], d, round(d / np.cos(np.radians((r[2] + r[3]) / 2)), 4))

    cgrid = {r[0]: grid(r, 0.03) for r in regions}
    wgrid = {r[0]: grid(r, 0.027) for r in regions if r[6] != "hrdps"}

    runs = {d: latest_hrrr(d, a.hours) for d in ("conus", "alaska")}
    print("HRRR runs", {d: f"{r:%Y-%m-%d %HZ}" for d, r in runs.items()})
    mappers = {}
    cloud = {r[0]: [] for r in regions}           # per region list of (t, low, mid, high)
    hwind = {r[0]: ([], [], [], [], []) for r in regions if r[6] != "hrdps"}
    start = min(runs.values())
    for h in range(0, a.hours + 1):
        t = start + dt.timedelta(hours=h)
        if t < now - dt.timedelta(hours=3):
            continue
        ts = int(t.timestamp())
        fields = {}
        for dom, run in runs.items():
            fh = int((t - run).total_seconds() // 3600)
            if 0 <= fh <= a.hours:
                try:
                    fields[dom] = hrrr_hour(dom, run, fh)
                except Exception as e:
                    print(f"  HRRR {dom} f{fh:02d}: {e}", file=sys.stderr)
        for r in regions:
            rid = r[0]
            g = cgrid[rid]
            layers = gfs_clouds_at(gfs, g, ts)          # coarse fill first, HRRR on top
            for dom in ("alaska", "conus"):
                if fields.get(dom) is None:
                    continue
                la, lo, f = fields[dom]
                key = (dom, rid, "c")
                if key not in mappers:
                    mappers[key] = Mapper(la, lo, g)
                m = mappers[key]
                for i, c in enumerate(("low", "mid", "high")):
                    v = m.apply(f[c])
                    layers[i] = np.where(np.isfinite(v), v, layers[i])
            cloud[rid].append((ts, layers))
            if r[6] in ("conus", "ak"):
                dom = "conus" if r[6] == "conus" else "alaska"
                if fields.get(dom) is None:
                    continue
                la, lo, f = fields[dom]
                key = (dom, rid, "w")
                if key not in mappers:
                    mappers[key] = Mapper(la, lo, wgrid[rid])
                m = mappers[key]
                T, U, V, G, P = hwind[rid]
                u, v = m.apply(f["u"]), m.apply(f["v"])
                T.append(ts); U.append(u); V.append(v)
                G.append(np.fmax(m.apply(f["gust"]), np.hypot(u, v)))
                P.append(np.maximum(m.apply(f["prate"]), 0) * 3600)
        if h % 6 == 0:
            print(f"  HRRR {t:%a %H}Z done", file=sys.stderr)

    gen = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    for r in regions:
        rid = r[0]
        raw = bytearray()
        for _, layers in cloud[rid]:
            for l in layers:
                raw += np.clip(np.round(np.nan_to_num(l)), 0, 100).astype(np.uint8).tobytes()
        header = dict(format="ttcloud", version=1, region=rid,
                      source="NOAA HRRR 3 km (CONUS and Alaska) with GFS 0.25° where HRRR doesn't reach",
                      license="US public domain (NOAA/NCEP)", generated=gen,
                      times=[t for t, _ in cloud[rid]], fields=["low", "mid", "high"], units="%", **cgrid[rid].meta())
        write_pack(os.path.join(a.out, f"{rid}.ttcloud"), b"HCLD", header, raw)
        if rid in hwind and hwind[rid][0]:
            T, U, V, G, P = hwind[rid]
            dom = "conus" if r[6] == "conus" else "alaska"
            write_wind_rain(r, a.out, wgrid[rid].meta(), T, U, V, G, P,
                            f"NOAA HRRR 3 km {dom}, run {runs[dom]:%Y-%m-%dT%HZ}", "US public domain (NOAA/NCEP)")
    print("clouds and HRRR wind written")

    # HRDPS regions in parallel (GeoMet WCS, four fields an hour each).
    s0, s1, ref = hrdps_range("HRDPS.CONTINENTAL_WSPD")
    g0, g1, _ = hrdps_range("HRDPS-WEonG_2.5km_WindGust")
    first, last = max(s0, g0), min(s1, g1)
    hstart = max(first, now - dt.timedelta(hours=3))
    hours = int((last - hstart).total_seconds() // 3600) + 1
    print(f"HRDPS run {ref}: {hours} hours")

    def one(r):
        with tempfile.TemporaryDirectory() as work:
            res = build_hrdps(r, hstart, hours, ref, work)
        write_wind_rain(r, a.out, *res)
        return r[0]

    with cf.ThreadPoolExecutor(4) as ex:
        for rid in ex.map(one, [r for r in regions if r[6] == "hrdps"]):
            print("  wrote", rid)

    manifest = dict(format="ttregions", version=1, generated=gen,
                    regions=[dict(id=r[0], name=r[1], south=r[2], north=r[3], west=r[4], east=r[5],
                                  wind=("HRDPS" if r[6] == "hrdps" else "HRRR"),
                                  files=[f"{r[0]}.ttwind", f"{r[0]}.ttrain", f"{r[0]}.ttcloud", f"{r[0]}.ttoutlook"])
                             for r in REGIONS])
    json.dump(manifest, open(os.path.join(a.out, "regions.json"), "w"), indent=1)
    for fn in sorted(os.listdir(a.out)):
        print(f"  {fn}: {os.path.getsize(os.path.join(a.out, fn)) / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
