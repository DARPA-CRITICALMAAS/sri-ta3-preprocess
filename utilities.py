import numpy as np
import rasterio
import pandas as pd
import geopandas as gpd
from sklearn import neighbors
from pathlib import Path
from tqdm import tqdm
import s2sphere as s2
from shapely.geometry import Polygon

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


def load_dataset(filename='data/2021_Table04_Datacube.csv', encoding_type='latin-1', index_col=None):
    df = pd.read_csv(filename, encoding=encoding_type, index_col=index_col)
    return df


def load_rasters(raster_files, rasters_path="data/LAWLEY22-RAW/geophysics/", verbosity=0):
    rasters = []
    for raster_file in raster_files:
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


def compute_bounds(rasters, region='Australia', option='-180_to_180', esp=1.0, verbosity=0):
    if option == '-180_to_180':
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
    elif option == 'compare_with_actual':
        if region == 'Australia':
            eyeballed_bounds = {'left':112.9, 'bottom':-43.6, 'right':153.6, 'top':-9.5} # including Tasmania but not Macquarie Island
            # eyeballed_bounds = {'left': 112.9, 'bottom': -54.5, 'right': 159, 'top': -9.5} # including Tasmania and Macquarie Island
        elif region == 'USCanada':
            eyeballed_bounds = {'left':-187.5, 'bottom':24.5, 'right':-52.6, 'top':83.15} # not including Hawaii
        bottom = left = 180.0
        top = right = -180.0
        # # iterating through rasters
        for raster in rasters:
            if raster.bounds.left < left and raster.bounds.left > (eyeballed_bounds['left'] - esp):
                left = raster.bounds.left
            if raster.bounds.bottom < bottom and raster.bounds.bottom > (eyeballed_bounds['bottom'] - esp): 
                bottom = raster.bounds.bottom
            if raster.bounds.top > top and raster.bounds.top < (eyeballed_bounds['top'] + esp):
                top = raster.bounds.top
            if raster.bounds.right > right and raster.bounds.right < (eyeballed_bounds['right'] + esp): 
                right = raster.bounds.right
    elif option == 'majority':
        # coming later
        raise NotImplementedError

    grid_bounds = {"bottom":bottom, "top":top, "left":left, "right":right}
    if verbosity:
        print(f"Computed bounds - {grid_bounds}\n")
    return grid_bounds


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


def s2_cell_center(cell_id):
    c1 = s2.Cell(cell_id)
    c0 = s2.LatLng.from_point(c1.get_center())
    return (c0.lng().degrees, c0.lat().degrees)


def s2_cell_polygon(cell_id):
    c1 = s2.Cell(cell_id)
    s2_vertices = [s2.LatLng.from_point(c1.get_vertex(i)) for i in range(4)]
    corr_vertices = [(s2_verx.lng().degrees, s2_verx.lat().degrees) for s2_verx in s2_vertices]
    corr_vertices.append((s2_vertices[0].lng().degrees, s2_vertices[0].lat().degrees))
    polygon = Polygon(corr_vertices)
    return polygon


def s2_cell_neighbors(c1, nn = 8):
    if nn == 4:
        neigh_id_list = [] # ???
    elif nn == 8:
        neigh_id_list = [neigh.id() for neigh in c1.get_all_neighbors(12)]
    return neigh_id_list


def init_datacube(initial_col, empty_cols, verbosity=0):
    datacube = pd.DataFrame(data=initial_col)
    for col in empty_cols:
        datacube[col] = np.nan
    datacube["mask"] = True
    if verbosity:
        datacube.head()
    return datacube


def mvt_dep_occur_to_s2cells(s2_df, mvt_df, colname='MVT_Deposit'):
    assert 's2_cell_id' in s2_df.columns.to_list()
    assert colname in ['MVT_Deposit', 'MVT_Occurrence']
    if colname not in s2_df.columns.to_list():
        s2_df[colname] = False
    list_s2_IDs = np.array(s2_df['s2_cell_id'])
    notrecog = []
    for index in range(len(mvt_df)):
        lat, lng = mvt_df.loc[index,'Latitude'], mvt_df.loc[index,'Longitude']
        s2_cell = latlong_to_cellid(lat, lng, s2_level=12)
        s2_cellid = s2_cell.id()
        if s2_cellid in list_s2_IDs:
            location = np.where(list_s2_IDs == s2_cellid)[0][0]
            s2_df.at[location, colname] = True
        else:
            notrecog.append(s2_cellid)
    return s2_df, notrecog


def process_raw_deposit_file(csv_file, csv_path='data/', region='Australia', dep_grp='MVT', verbosity=0):
    assert region in ['Australia', 'USCanada']
    assert dep_grp in ['MVT', 'CD']
    df = load_dataset(f'{csv_path}{csv_file}')

    if dep_grp == 'MVT':
        df = df[df['Dep_Grp'] == 'Mississippi Valley-type (Zn, Pb)']
    elif dep_grp == 'CD':
        df = df[df['Dep_Grp'] == 'Clastic-dominated (Zn, Pb)']
    df = df[df['Site_Class'] != 'District'] # only Deposit and Occurrence

    if region == 'Australia':
        df = df[df['Admin'] == 'Australia']
    elif region == 'USCanada':
        df = df[df['Admin'] != 'Australia']

    df = df.drop(columns=['Admin','Dep_Grp','Source','Dep_Type','Dep_Name','Tonnage_Mt','Cu_pct','Zn_pct','Pb_pct','Ag_ppm','Au_ppm'])

    df_deposit = df[df['Site_Class']=='Deposit']
    df_deposit = df_deposit.drop_duplicates().reset_index(drop=True)
    df_occurrence = df[df['Site_Class']=='Occurrence']
    df_occurrence = df_occurrence.drop_duplicates().reset_index(drop=True)
    return df_deposit, df_occurrence


def neighbor_deposits(df, deptype='MVT'):
    assert deptype in ['MVT','CD']
    s2_cells_all = np.array(df['s2_cell_id'])

    df[f'{deptype}_DepositOccurrence'] = df.apply(lambda row: True if True in [row[f'{deptype}_Deposit'], row[f'{deptype}_Occurrence']] else False, axis=1)
    s2_cells_do = np.array(df[df[f'{deptype}_DepositOccurrence']==True]['s2_cell_id'])

    s2_cells_don = s2_cells_do
    for cell in s2_cells_do:
        neighbors_cells = s2_cell_neighbors(s2id_to_cellid(int(cell)), nn=8)
        s2_cells_don = np.append(s2_cells_don, np.array(neighbors_cells))
        # s2_cells_all = np.concatenate([s2_cells_all, neighbors_cells])

    s2_cells_don = np.unique(s2_cells_don)
    s2_cells_don = np.intersect1d(s2_cells_don, s2_cells_all)
    df[f'{deptype}_DepositOccurrenceNeighbors'] = False
    for cell in s2_cells_don:
        location = np.where(s2_cells_all == cell)[0][0]
        df.at[location, f'{deptype}_DepositOccurrenceNeighbors'] = True
    return df