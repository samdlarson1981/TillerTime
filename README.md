# Tiller Time: forecast packs

This repo builds the current packs for **Tiller Time**, an iPhone chartplotter for the Salish Sea.

Four times a day a GitHub Actions job builds two small forecast files and publishes them here:

| Layer | File | Source | Grid / range |
|---|---|---|---|
| Currents | [salishsea.ttcur](https://github.com/samdlarson1981/TillerTime/releases/download/currents-latest/salishsea.ttcur) | SalishSeaCast (UBC) | ~450 m, hourly, ~35–40 h ahead, ~5 MB |
| Wind and gusts | [salishsea.ttwind](https://github.com/samdlarson1981/TillerTime/releases/download/currents-latest/salishsea.ttwind) | ECCC HRDPS via MSC GeoMet | 2.5 km, hourly, ~45 h ahead, ~3 MB |

Tiller Time downloads both automatically while it's open. You can also tap **More › Forecasts › Update now**.

See [tools/currents/README.md](tools/currents/README.md) for how the pack is built, the file format and known limits.

**Credits:** Wind contains information licensed under the Open Government Licence – Canada (Environment and Climate Change Canada, HRDPS). SalishSeaCast NEMO model results are © the SalishSeaCast Project Contributors and The University of British Columbia, licensed under Apache 2.0. Cite Soontiens et al. (2016) and Soontiens and Allen (2017).

**Not for navigation on its own.** These are model forecasts. They underestimate the flow in passes narrower than the grid (Deception Pass, Active Pass, Dodd Narrows and others).
