SPATIAL LAYERS
==============
The two Kuskokwim spatial layers required by the mapping, clustering, and
cluster-analysis stages are stored in the "Spatial Data" subfolder (copied in
from the original analysis). Each shapefile includes its sidecar files
(.shp, .shx, .dbf, .prj, ...):

    data/spatial/Spatial Data/KuskoUSGS_HUC_joined.shp     (stream edges)
    data/spatial/Spatial Data/Kusko_basin.shp              (basin outline)

These originate from the French et al. (2020) Kuskokwim strontium isoscape
(see ../../README.md). The paths are set in ../../config.R (KUSKO_EDGES /
KUSKO_BASIN); edit those if you ever move the files.
