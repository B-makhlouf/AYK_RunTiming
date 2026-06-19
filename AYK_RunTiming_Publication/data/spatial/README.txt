SPATIAL LAYERS — KEEP THE ORIGINAL FOLDER STRUCTURE
===================================================
The two Kuskokwim spatial layers are required by the mapping, clustering, and
cluster-analysis stages. They are not redistributed in this repository. Copy
them in using the SAME structure as the original analysis (each shapefile needs
all of its sidecar files: .shp, .shx, .dbf, .prj, ...):

    data/spatial/Spatial Data/KuskoUSGS_HUC_joined.shp     (stream edges)
    data/spatial/isoscapes_new/Kusko/Kusko_basin.shp       (basin outline)

These mirror the original locations:
    .../Spatial Data/KuskoUSGS_HUC_joined.shp
    .../Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp

so you can drop the original "Spatial Data" and "isoscapes_new" folders in here
without renaming anything. If your paths differ, edit KUSKO_EDGES / KUSKO_BASIN
in ../../config.R instead.

Source: French et al. (2020) Kuskokwim strontium isoscape (see ../../README.md).
