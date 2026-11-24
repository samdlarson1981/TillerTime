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
