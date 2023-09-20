## 1. Implement preprocessing needed for woe_baseline (i.e. geophysical input variables)
- Functions needed: load raster data (`.tif`), load vector data (`.shp`), converting vector data (`.shp`) to raster data (`.tif`), "unifying" raster data (`.tif`),
- Existig tools to implement above functions:

| Tool | Supported Functions                                           |
|------|---------------------------------------------------------------|
| [rasterio](https://rasterio.readthedocs.io/en/stable/) | load raster data |
| [geopandas](https://geopandas.org/en/stable/getting_started/introduction.html#) | load vector data |
| [eis_toolkit](https://github.com/GispoCoding/eis_toolkit/tree/master) | converting vector data to raster data, "unifying" raster data |
| [s2sphere](https://s2sphere.readthedocs.io/en/latest/index.html) | Python implementation of S2 geometry library |
- Subtasks:
    - Make repo that submodules the eis_toolkit, and install rasterio and geopandas (done)
    - Make notebook to *"Implement preprocessing needed for geophysical input variables"* (done)
    - Write notebook cell(s) that load geophysical raster files (`.tif`). (done)
    - Write notebook cell(s) that load geophysical vector files (`.shp`). Includes grav and mag worms. (done)
    - Write notebook cell(s) that convert geophysical vector files into appropriate raster files. E.g. grav worm vector files need to be converted to proximity.
    - Write notebook cell(s) that convert training vector files into raster files. E.g. MVT mineral occurences with Lat/Long into a mask raster (presence/absence).
    - Write notebook cell(s) that unifies all geophysical and training raster data
    

1. Gather / load all vector (.shp), raster (.tif), and deposit/occurences (.csv) data relevant to the experiment (done)
2. Define bounds - choose the max size (left-most, right-most, top-most, bottom-most) bound for all vector and raster files (done)
3. Extract the grid - query the S2 library for the grid cells at appropriate resolution within the chosen bounds (Australia - done, USCanada - done)
4. Now initialize the "datacube" - make a geopandas [or pandas?] table where each row is one cell with the grid from 3. and has no columns yet (Australia - done, USCanada - done)
5. Add grid relevant info to data cube - make first column(s) be address (and other relevant, reusable info, e.g. lat/log) of every grid cells (Australia - done, USCanada - done | NOTE: only kept S2_ID name, any other infocan be easily recovered with s2sphere)
6. Add mask to datacube - make next column be for mask of whether data is present or not
7. Add all other columns to the table with empty values - i.e. the input variables from .shp or .tif
8. Populate datacube with info
   a. Start a for loop where one loop goes through grid cells and other loop goes through input data (columns of datacube which are vector / raster files above)
   b. for every grid cell, grab the boundaries of the grid cell
   c. --- for every input data, query that file (tif or shp) for the data value(s) within the grid cell boundaries from b. - there is likely a QGIS tool to do this
   d. once data is written to a grid cell, unmask it
9. output the datacube as raster or geopandas table
