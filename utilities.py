import numpy as np
import rasterio
import geopandas as gpd
from sklearn import neighbors
from pathlib import Path
from tqdm import tqdm
import s2sphere as s2

def get_input_var_files(region):
    if "aus" in region.lower():
        tif_files = [
            "GeophysicsLAB_Australia",
            "GeophysicsMoho_Australia",
            "GeophysicsSatelliteGravity_ShapeIndex_Australia",
            "GeophysicsGravity_Australia",
            "GeophysicsGravity_HGM_Australia",
            "GeophysicsGravity_UpCont30km_Australia",
            "GeophysicsGravity_UpCont30km_HGM_Australia",
            "GeophysicsMagRTP_Australia",
            "GeophysicsMagRTP_DeepSources_Australia",
            "GeophysicsMagRTP_VD_Australia",
            "GeophysicsMagRTP_HGM_Australia",
            "GeophysicsMagRTP_HGMDeepSources_Australia",
        ]
        
        shp_files = [
            "ShallowGravitySources_Worms_Australia",
            "DeepGravitySources_Worms_Australia",
            "ShallowMagSources_Worms_Australia",
            "DeepMagSources_Worms_Australia",
        ]
    else:
        tif_files = [
            "GeophysicsLAB_USCanada",
            "USCanada_Moho",
            "GeophysicsSatelliteGravity_ShapeIndex_USCanada",
            "GeophysicsGravity_USCanada",
            "GeophysicsGravity_HGM_USCanada",
            "GeophysicsGravity_Up30km_USCanada",
            "GeophysicsGravity_Up30km_HGM_USCanada",
            "GeophysicsMag_RTP_USCanada",
            "USCanadaMagRTP_DeepSources",
            "GeophysicsMag_RTP_VD_USCanada",
            "GeophysicsMag_RTP_HGM_USCanada",
            "USCanadaMagRTP_HGMDeepSources",
        ]

        shp_files = [
            "ShallowGravitySources_Worms_USCanada",
            "DeepGravitySources_Worms_USCanada",
            "ShallowMagSources_Worms_USCanada",
            "DeepMagSources_Worms_USCanada",
        ]
    return tif_files, shp_files


def load_rasters(raster_files, rasters_path="data/LAWLEY22-RAW/geophysics/", verbosity=0):
    rasters = []
    pbar = tqdm(raster_files)
    for raster_file in pbar:
        pbar.set_description(f"Loading {raster_file}")
        rasters.append(load_raster(raster_file, rasters_path, verbosity))
    return rasters


def load_raster(raster_file, raster_path="data/LAWLEY22-RAW/geophysics/", verbosity=0):
    raster = rasterio.open(f"{raster_path}{raster_file}.tif")
    if verbosity:
        print(f"-------- {raster_file} raster details --------\n")
        info = {i: dtype for i, dtype in zip(raster.indexes, raster.dtypes)}
        print(f"Raster bands and dtypes:\n{info}\n\n")
        print(f"Coordinate reference system:\n{raster.crs}\n\n")
        print(f"Bounds:{raster.bounds},Size:{raster.shape},Resolution:{raster.res}\n\n")
    return raster


def load_vectors(vector_files, vectors_path="data/LAWLEY22-RAW/geophysics/", verbosity=0):
    vectors = []
    pbar = tqdm(vector_files)
    for vector_file in pbar:
        pbar.set_description(f"Loading {vector_file}")
        vectors.append(load_vector(vector_file, vectors_path, verbosity))
    return vectors


def load_vector(vector_file, vector_path="data/LAWLEY22-RAW/geophysics/", verbosity=0):
    vector = gpd.read_file(f"{vector_path}{vector_file}.shp")
    if verbosity:
        print(f"-------- {vector_file} raster details --------\n")
        print(f"{vector.head()}\n\n")
    return vector

def get_vector_pts(vector):
    try:
        vector_pts = vector.loc[:,["Longitude", "Latitude"]].values
    except:
        vector_pts = vector.loc[:,["Lon_WGS", "Lat_WGS"]].values
    return vector_pts

def proximity_raster_of_vector_points(raster, vector_file, vector):
    # get the lat / lons of vector points
    vector_pts = get_vector_pts(vector)
    # extract the lat / lons of raster
    rows, cols = np.mgrid[0:raster.height:1,0:raster.width:1].reshape((-1, raster.width*raster.height))
    raster_pts = np.dot(np.asarray(raster.transform.column_vectors).T, np.vstack((cols, rows, np.ones_like(rows)))).T
    # call neighbors computation - slowest part of code, O(NlogN)
    tree = neighbors.KDTree(vector_pts, metric="euclidean")
    proximity, _ = tree.query(raster_pts, k=1)
    proximity = np.asarray(proximity).reshape(raster.height, raster.width)
    proximity = np.ma.array(proximity, mask=~raster.read_masks(1).astype(bool)).filled(raster.nodata)
    # clone above raster with same mask, updated values to proximity
    vector_tif_file = str(Path(raster.files[0]).parent / Path(vector_file + ".tif"))
    vector_raster = rasterio.open(vector_tif_file, "w", **raster.meta)
    vector_raster.write(proximity, 1)
    vector_raster.close()


def compute_bounds(rasters):
    bottom = left = 180.0
    top = right = -180.0
    for raster in rasters:
        if raster.bounds.left < left:
            left = raster.bounds.left
        if raster.bounds.bottom < bottom:
            bottom = raster.bounds.bottom
        if raster.bounds.top > top:
            top = raster.bounds.top
        if raster.bounds.right > right:
            right = raster.bounds.right
    return {"bottom":bottom, "top":top, "left":left, "right":right}


def s2id_to_cellid(id):
    # input: id - must be an integer
    # ouput: S2 CellId, where level is automatically determined
    return s2.CellId(id)


def latlong_to_cellid(lat, lon, s2_level):
    # takes (lat,lon) point, and finds the CellId w/ respect to s2_level
    point = s2.LatLng.from_degrees(lat, lon)
    cellid = s2.CellId.from_lat_lng(point).parent(s2_level)
    return cellid


def region_of_cellids(dict_bounds, s2_level):
    # takes bounds_australia (WORKS) or bounds_uscanada (ISSUES)
    # outputs a list of all CellIds that cover the region
    coverer = s2.RegionCoverer()
    coverer.min_level = coverer.max_level = s2_level

    # order matters, particularly for the s2.LatLngRect() function
    min_lat_lng = s2.LatLng.from_degrees(dict_bounds['bottom'], dict_bounds['left']) # bottom left corner
    max_lat_lng = s2.LatLng.from_degrees(dict_bounds['top'], dict_bounds['right']) # top right corner

    # Create an S2LatLngRect from the boundary points
    rect = s2.LatLngRect(min_lat_lng, max_lat_lng)

    cell_ids = coverer.get_covering(rect)
    return cell_ids
