# Tiller Time: forecast packs

This repo builds the current packs for **Tiller Time**, an iPhone chartplotter for the Salish Sea.

Four times a day a GitHub Actions job builds small forecast files and publishes them here:

| Layer | File | Source | Grid / range |
|---|---|---|---|
| Currents | [salishsea.ttcur](https://github.com/samdlarson1981/TillerTime/releases/download/currents-latest/salishsea.ttcur) | SalishSeaCast (UBC) | ~450 m, hourly, ~35–40 h ahead, ~5 MB |
| Wind and gusts | [salishsea.ttwind](https://github.com/samdlarson1981/TillerTime/releases/download/currents-latest/salishsea.ttwind) | ECCC HRDPS via MSC GeoMet | 2.5 km, hourly, ~45 h ahead, ~3 MB |
| Rain | [salishsea.ttrain](https://github.com/samdlarson1981/TillerTime/releases/download/currents-latest/salishsea.ttrain) | ECCC HRDPS precipitation rate | 2.5 km, hourly, ~45 h ahead |
| Clouds (low / mid / high) | [salishsea.ttcloud](https://github.com/samdlarson1981/TillerTime/releases/download/currents-latest/salishsea.ttcloud) | NOAA HRRR | 3 km, hourly, ~48 h |
| 10-day outlook | [salishsea.ttoutlook](https://github.com/samdlarson1981/TillerTime/releases/download/currents-latest/salishsea.ttoutlook) | NOAA GFS | 0.25°, 3-hourly, 10 days: clouds, wind, gusts, rain |
| Tides | [salishsea.tides.json](https://github.com/samdlarson1981/TillerTime/releases/download/currents-latest/salishsea.tides.json) | NOAA CO-OPS + Canadian Hydrographic Service | ~305 stations, highs and lows for 8 days |

Tiller Time downloads them automatically while it's open. You can also tap **More › Forecasts › Update now**.

See [tools/currents/README.md](tools/currents/README.md) for how the pack is built, the file format and known limits.

**Credits:** Tide predictions: NOAA CO-OPS and the Canadian Hydrographic Service. Wind and rain contain information licensed under the Open Government Licence – Canada (Environment and Climate Change Canada, HRDPS). SalishSeaCast NEMO model results are © the SalishSeaCast Project Contributors and The University of British Columbia, licensed under Apache 2.0. Cite Soontiens et al. (2016) and Soontiens and Allen (2017).

**Not for navigation on its own.** These are model forecasts. They underestimate the flow in passes narrower than the grid (Deception Pass, Active Pass, Dodd Narrows and others).

## Vector charts

`tools/enc/build_enc_pack.py` turns NOAA's free electronic navigational charts (S-57 ENCs,
public domain) for the US Salish Sea into one compact pack, `charts_salish.ttenc`, that the
app draws on the phone: depth areas, contours, soundings, land, docks, rocks, wrecks,
buoys, beacons, lights and place names. `.github/workflows/charts.yml` rebuilds it every
Monday. Canadian waters use the CHS chart service in the app. Not for navigation.

### BC vector charts

`tools/enc/build_bc_pack.py` makes a free chart for BC waters (Victoria to Desolation Sound's
southern approaches) from CHS NONNA-10 survey depths (Open Government Licence – Canada) and
OpenStreetMap land, buoys, lights, docks and names (© OpenStreetMap contributors, ODbL).
CHS's data service refuses GitHub's cloud runners, so the depth stage runs elsewhere and is
committed as `data/bc_depths.henc`; `.github/workflows/charts-bc.yml` adds the OpenStreetMap
layers monthly and publishes `charts_bc.ttenc`. Not for navigation.
