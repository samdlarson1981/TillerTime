# Tiller Time: current forecasts

This repo builds the current packs for **Tiller Time**, an iPhone chartplotter for the Salish Sea.

Four times a day a GitHub Actions job takes the latest **SalishSeaCast** forecast (University of British Columbia) and turns it into a ~5 MB `salishsea.ttcur` file: near-surface currents on a ~450 m grid, hourly, for about 35–40 hours ahead. The job publishes it here:

https://github.com/samdlarson1981/TillerTime/releases/download/currents-latest/salishsea.ttcur

In Tiller Time, paste that address into **More › Currents › Pack address** and tap **Download latest**.

See [tools/currents/README.md](tools/currents/README.md) for how the pack is built, the file format and known limits.

**Credits:** SalishSeaCast NEMO model results are © the SalishSeaCast Project Contributors and The University of British Columbia, licensed under Apache 2.0. Cite Soontiens et al. (2016) and Soontiens and Allen (2017).

**Not for navigation on its own.** These are model forecasts. They underestimate the flow in passes narrower than the grid (Deception Pass, Active Pass, Dodd Narrows and others).
