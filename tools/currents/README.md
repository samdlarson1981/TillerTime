# Current packs

`build_current_pack.py` turns the latest **SalishSeaCast** forecast (University of British Columbia) into a small `.ttcur` file that Tiller Time draws as current arrows.

- **Coverage:** Juan de Fuca Strait, Strait of Georgia, San Juan and Gulf Islands, Puget Sound, up to Johnstone Strait.
- **Grid:** about 440 × 500 m, near-surface (upper 5 m average), hourly.
- **Time span:** from 3 hours ago to the end of the latest run, about 35–40 hours.
- **Size:** about 5 MB.
- **Licence:** the model results are © SalishSeaCast Project Contributors and UBC, licensed under Apache 2.0. Credit them when sharing.

## Make a pack on your Mac

```
pip3 install numpy scipy
python3 tools/currents/build_current_pack.py -o salishsea.ttcur
```

AirDrop `salishsea.ttcur` to your iPhone and choose **Tiller Time**. You can also open it from Files, or use Tiller Time › More › Currents › Import.

## Automatic packs (free)

`.github/workflows/currents.yml` rebuilds the pack four times a day on GitHub Actions and publishes it at:

`https://github.com/<you>/<repo>/releases/download/currents-latest/salishsea.ttcur`

Put that address in Tiller Time › More › Currents › Pack address, then tap **Download latest** before you leave the dock.

## Known limits

- The grid is coarser than narrow passes: Deception Pass, Active Pass, Porlier, Gabriola, Dodd Narrows, the Yucultas, Seymour Narrows and others. Peak flow there is badly underestimated, and Tiller Time shows a warning near them. Use the tide tables for those passes.
- The values are a model forecast, not observations. Compare them with what you see on the water.

## File format (version 1)

```
"HCUR"  u16 version=1  u16 reserved  u32 headerLength
header JSON: lat0 lon0 dlat dlon nlat nlon nWater times[] scale rawBytes source ...
raw DEFLATE block (no zlib header), inflating to rawBytes:
    int32 index[nlat*nlon]      regular cell -> water cell, -1 = land   (row = latitude)
    for each time: int16 u[nWater], int16 v[nWater]                     (× scale = m/s, east/north)
```

The model grid is curvilinear, so it is interpolated onto the regular grid using triangles, and no triangle may span more than 1.2 km. That keeps currents from leaking across islands and spits. Checked against the model output, the median error is 0.012 m/s.
