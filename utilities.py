import rasterio
from rasterio.transform import from_origin
import geopandas as gpd
import numpy as np
import pandas as pd
import s2sphere as s2

# computed bounds for Australia and US/Canada (Note: seems like US/Canada bounds have issues)
bounds_australia = {'left': 112.91657, 'bottom': -54.7421128072652, 'right': 159.100865, 'top': -9.221274999999999}
bounds_uscanada = {'left': -191.00025, 'bottom': 15.696449999999999, 'right': 179.754051886584, 'top': 90.05655}

def load_dataset(filename='data/2021_Table04_Datacube.csv', encoding_type='latin-1', index_col=None):
    print(f"Loading {filename}")
    df = pd.read_csv(filename, encoding=encoding_type, index_col=index_col)
    return df

def load_tif_file(tif_filename):
    # requires rasterio library
    print(f"Loading {tif_filename}")
    return rasterio.open(tif_filename)

def load_shp_file(shp_filename):
    # requires rasterio library
    print(f"Loading {shp_filename}")
    return gpd.read_file(shp_filename)

def determine_parameters(df, kmpdeg_lat0=111.1, crs = 'EPSG:4326', kminpix=2.56318,
                                spec_resolution=None, spec_latitude=None, ):
    # function to determine parameteres for futher tif file construction

    #### --------- INPUTS ---------
    # - kmpdeg_lat0: kilometers per degree on Equator, default: 111.1
    # - crs: coordinate reference system, default: EPSG:4326
    # - kminpix: kilometers in pixel in outcome/produced grid, default: 2.56318

    # - spec_resolution: specified resolution, default: None
    # - spec_latitude: specified latitude, default: None

    ### --------- OUTPUT(S) ---------
    # - dictionary with the following parameters: 'latitude_name' 'longitude_name' 'resolution' 
    #                                               'resolution_latitude' 'height' 'width' 'boundaries'
    #                                               'transform' 'np_transform' 'no_data_constant'
    
    no_data_constant = -340282346638528859811704183484516925440.0

    columns = df.columns.to_list()
    latitude_col_name = [col for col in columns if 'lat' in col.lower()][0] # e.g., 'Lat_WGS84' or 'Lat_WGS'
    longitude_col_name = [col for col in columns if 'lon' in col.lower()][0] # e.g., 'Lon_WGS84' or 'Lon_WGS'

    bounds = {
        'left' : df[longitude_col_name].min(), 
        'bottom' : df[latitude_col_name].min(), 
        'right' : df[longitude_col_name].max(), 
        'top' : df[latitude_col_name].max()
    }

    # define resolution
    if spec_resolution:
        resolution = spec_resolution
    else:
        if spec_latitude:
            latitude = spec_latitude
        else:
            # we can use the average latitude of the region
            latitude = (bounds['top'] + bounds['bottom']) / 2.0
    
        km_in_a_degree = kmpdeg_lat0 * np.cos(np.deg2rad(latitude))
        resolution = kminpix / km_in_a_degree

    height = int((bounds['top'] - bounds['bottom']) / resolution)
    width = int((bounds['right'] - bounds['left']) / resolution)

    transform = from_origin(bounds['left']-resolution/2, bounds['top']-resolution/2, resolution, resolution)

    np_matrix_transform = np.array(transform).reshape(3,3)

    output_dict = {
        'latitude_name' : latitude_col_name,
        'longitude_name' : longitude_col_name,
        'resolution' : resolution,
        'resolution_latitude' : latitude,
        'height' : height,
        'width' : width,
        'boundaries' : bounds,
        'transform' : transform,
        'np_transform' : np_matrix_transform,
        'no_data_constant' : no_data_constant,
        'crs' : crs,
    }
    return output_dict

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
