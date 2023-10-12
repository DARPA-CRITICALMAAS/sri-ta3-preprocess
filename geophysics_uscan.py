# %% [markdown]
# Imports and set relevant paths

# %%
import numpy as np
from tqdm import tqdm
import rasterio

import utilities as utils

RAW_DATA_DIR = "data/LAWLEY22-RAW/"
DERIV_DATA_DIR = "data/LAWLEY22-DERIV/"

# %% [markdown]
# Load the names of ScienceBase raster (`.tif`) and vector (`.shp`) files.

# %%
tifs, shps = utils.get_input_var_files("USCAN")

# %% [markdown]
# Loads the raster data

# %%
rasters = utils.load_rasters(tifs, rasters_path=RAW_DATA_DIR, verbosity=1)

# %% [markdown]
# We'll need to upsample rasters that have too low a resolution

# %%
target_resolution = (0.008, 0.008)
rasters = utils.resample_rasters(target_resolution, rasters, tifs)

# %% [markdown]
# Loads rasters of the vector data if available; otherwise generates them

# %%
try:
    rasters += utils.load_rasters(shps, rasters_path=RAW_DATA_DIR, verbosity=1)
except rasterio.RasterioIOError:
    base_raster = rasters[-1] # defaults to intermediate resolution raster
    vectors = utils.load_vectors(shps, vectors_path=RAW_DATA_DIR, verbosity=0)
    pbar = tqdm(zip(shps, vectors))
    for shp, vector in pbar:
        pbar.set_description(f"Processing {shp}")
        utils.proximity_raster_of_vector_points(base_raster, shp, vector)
    rasters += utils.load_rasters(shps, rasters_path=RAW_DATA_DIR, verbosity=1)

# %% [markdown]
# Loads the base grid for all data if available; otherwise generates it

# %%
grid_cell_ids = utils.generate_s2_grid(rasters, DERIV_DATA_DIR, "s2_grid_uscan")

# %% [markdown]
# Initialize the datacube

# %%
datacube = utils.init_datacube({"s2_cell_id": grid_cell_ids}, ["s2_cell_center", "s2_cell_poly"] + tifs + shps, verbosity=1)

# %% [markdown]
# Load file with Deposits and Occurrences

# %%
df_dep, df_occ = utils.process_raw_deposit_file('GeologyMineralOccurrences_USCanada_Australia.csv', csv_path='data/LAWLEY22-RAW/', region='USCanada', dep_grp='MVT')

# %% [markdown]
# Adding MVT_Deposit, MVT_Occurrence columns to datacube

# %%
datacube, notrecogdep = utils.mvt_dep_occur_to_s2cells(datacube, df_dep, colname='MVT_Deposit')
datacube, notrecogocc = utils.mvt_dep_occur_to_s2cells(datacube, df_occ, colname='MVT_Occurrence')

# %% [markdown]
# Add neighbors column to datacube

# %%
datacube = utils.neighbor_deposits(datacube, deptype='MVT')

# %% [markdown]
# Populate the datacube using as many process as available CPUs

# %%
datacube = utils.populate_datacube(datacube, tifs+shps)

# %% [markdown]
# Final filtering of s2 cells that contain no data

# %%
print(f"Removing {np.count_nonzero(np.isnan(datacube.loc[:, tifs+shps].values).all(axis=1))} s2 cells that have no geophysical data")
datacube = datacube[~np.isnan(datacube.loc[:, tifs+shps].values).all(axis=1)]

# %% [markdown]
# Store the datacube for future use

# %%
datacube.to_csv(f"{DERIV_DATA_DIR}datacube_uscan.csv", index=False)


