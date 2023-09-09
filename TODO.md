## 1. Implement preprocessing needed for woe_baseline (i.e. geophysical input variables)
- Functions needed: load raster data (`.tif`), load vector data (`.shp`), converting vector data (`.shp`) to raster data (`.tif`), "unifying" raster data (`.tif`),
- Existig tools to implement above functions:

| Tool | Supported Functions                                           |
|------|---------------------------------------------------------------|
| [rasterio](https://rasterio.readthedocs.io/en/stable/) | load raster data |
| [geopandas](https://geopandas.org/en/stable/getting_started/introduction.html#) | load vector data |
| [eis_toolkit](https://github.com/GispoCoding/eis_toolkit/tree/master) | converting vector data to raster data, "unifying" raster data |
- Subtasks:
    - Make repo that submodules the eis_toolkit, and install rasterio and geopandas (done)
    - Make notebook to *"Implement preprocessing needed for geophysical input variables"* (done)
    - Write notebook cell(s) that load geophysical raster files (`.tif`). (done)
    - Write notebook cell(s) that load geophysical vector files (`.shp`). Includes grav and mag worms. (done)
    - Write notebook cell(s) that convert geophysical vector files into appropriate raster files. E.g. grav worm vector files need to be converted to proximity.
    - Write notebook cell(s) that convert training vector files into raster files. E.g. MVT mineral occurences with Lat/Long into a mask raster (presence/absence).
    - Write notebook cell(s) that unifies all geophysical and training raster data
    