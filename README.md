### Install
1. Create the conda environment - `conda create -n ta3-preprocessing python=3.8`
2. Activate conda environment - `conda activate ta3-preprocessing`
3. Install gdal and raterio using conda - `conda install -c conda-forge "gdal>=3.0.2, <4.0.0" rasterio`
3. Install tool libraries - `pip3 install -r requirements.txt`


### Prepare input datasets
1. Download, extract, and place all relevant datasets ("Child Items" {1, [7-13]} as published on [the Science base page](https://www.sciencebase.gov/catalog/item/6193e9f3d34eb622f68f13a5)) from [Lawley'22](https://www.sciencedirect.com/science/article/pii/S0169136821006612) into the `data/LAWLEY22-RAW` folder following the directory structure below .
2. Download, extract, and place the "Ocean" [vector](https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_ocean.zip) from the [highest resolution Natural Earth](https://www.naturalearthdata.com/downloads/10m-physical-vectors/) free vectors into the `data/OCEAN` folder following the directory strcture below.

```
data
├── LAWLEY22-RAW # contains all relevant datasets from ScienceBase
│   ├── GeophysicsGravity_Australia.tif 
│   ├── GeophysicsGravity_HGM_Australia.tif
│   └── ...
├── LAWLEY22-DERIV # will contain all generated data
│   ├── s2_grid_aus.npy # provided for convenience
│   ├── s2_grid_uscan.npy  # provided for convenience
│   └── ...
└── OCEAN # contains the ocean vector data
    ├── ne_10m_ocean.shp
    └── ...
```

### Preprocessing Geophysical Data
Preprocessing to convert geophysical and Basin-hosted (CD/SEDEX and MVT) Zn-Pb deposits/prospects datasets into a datacube using the [S2 grid system](https://s2geometry.io/) for all continents in Lawley'22 are within a corresponding notebook. Please open and view/run the cells of the notebook to view/run the preprocessing steps. For convenience, the `generate_datacubes.sh` script is also provided to generate datacubes for all continents using the s2 cell grid system.

