import rasterio
from rasterio.transform import from_origin
import geopandas as gpd
import numpy as np
import pandas as pd

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
    #                                               'transform' 'np_transform' 'no_data_constant' 'crs'
    
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
