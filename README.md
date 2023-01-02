# Data cube for Haskell University Workshop

![ESIIL_datacube](https://user-images.githubusercontent.com/67020853/210274189-ec816a81-b69b-41dd-a731-b7f163fa4de0.png)

## 1. Create Area of Interest (Python code)

```python
# Import packages 
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon
import fiona
from fiona.crs import from_epsg

# Create an empty geopandas GeoDataFrame
newdata1 = gpd.GeoDataFrame()
newdata1

# Create field geometry
newdata1['geometry'] = None
newdata1

# Create a Shapely polygon from the coordinate
coordinates1 = [(-95.16662, 38.90447), (-95.16662, 39.03350), (-95.34454, 39.03350), (-95.34454, 38.90447)]

extent1 = Polygon(coordinates1)
extent1 

# Insert the polygon into 'geometry' -column at index 0
newdata1.loc[0, 'geometry'] = extent1

# Add a new column and insert data
newdata1.loc[0, 'Location'] = 'Haskell_AOI'

#Add crs (WGS84)
newdata1.crs = from_epsg(4326)
newdata1.crs

# Determine the output path for the Shapefile
outfp1 = r"..\Haskell_wgs84.shp"

# Determine the output path for the GeoJSON
outfp2 = r"..\Haskell_wgs84.geojson"

# Write the data into that Shapefile
newdata1.to_file(outfp1)

# Write the data into GeoJSON
newdata1.to_file(outfp2, driver="GeoJSON")
```

## 2. Download data

#### NLCD: National Land Cover Database 2019, Landcover & Imperviousness (NLCD2019)

Download from https://www.mrlc.gov/data?f%5B0%5D=region%3Aconus&f%5B1%5D=year%3A2019

#### SRTM: NASA Shuttle Radar Topography Mission Global 1 arc second V003

Download from https://search.earthdata.nasa.gov/search

*Use the geojson file generated in section 1 for search

#### ECOSTRESS: ECOSTRESS Water Use Efficiency Daily L4 Global 70m V001 & ECOSTRESS Geolocation Daily L1B Global 70m V001

Download from https://search.earthdata.nasa.gov/search?q=ECOSTRESS

*Use the shape file generated in section 1 for search

**Restrict the L4 data search by year, season, and hours of the day

***Search and download L1B data by using the relative L4 data orbits 

#### GEDI: GEDI L1B Geolocated Waveform Data Global Footprint Level V002 & GEDI L2A Elevation and Height Metrics Data Global Footprint Level V002

Download using the rGEDI package (R code): https://github.com/carlos-alberto-silva/rGEDI
```r
# Set AOI box coordinates
ul_lat<- 39.03350
lr_lat<- 38.90447
ul_lon<- -95.34454
lr_lon<- -95.16662

# Specify the date range
daterange=c("2020-01-01","2020-12-31")

# Get path to GEDI data
gLevel1B<-gedifinder(product="GEDI01_B",ul_lat, ul_lon, lr_lat, lr_lon,version="002",daterange=daterange)
gLevel2A<-gedifinder(product="GEDI02_A",ul_lat, ul_lon, lr_lat, lr_lon,version="002",daterange=daterange)
```
## 3. Mosaic scenes (SRTM) (Python code)

```python=
# Import packages
import os
from osgeo import gdal
import numpy as np

# List of .htg filenames
filenames = ['N38W096.htg', 'N39W096.htg']

# Read in the .htg files
datasets = [gdal.Open(f) for f in filenames]

# Create a mosaic
mosaic_SRTM, _ = gdal.Composite(output, datasets)

# Save the mosaic as a .tiff file
driver = gdal.GetDriverByName('GTiff')
mosaic.CreateCopy('mosaic_SRTM.tiff', driver)
```
### 4. Transform swath to grid (ECOSTRESS) (Python code)

Use script available at NASA EarthData Bitbucket: https://git.earthdata.nasa.gov/projects/LPDUR/repos/ecostress_swath2grid/browse

After setting-up the conda environment, run it in terminal as follows:

```python
> python ECOSTRESS_swath2grid.py --proj GEO --dir C:\Users\ECOSTRESS\
```
*If running the script in Windows/Git bash, type --dir path using forward slash

## 5. Create composite image (ECOSTRESS WUE) (Python code)

```python
import gdal
import numpy as np
import glob

# Get a list of all the geotif images in the specified directory
image_list = glob.glob('path/to/images/*_WUEavg_GEO.tif')

# Read the data from each image and store it in a list
image_data = []
for image in image_list:
    dataset = gdal.Open(image)
    image_data.append(dataset.ReadAsArray())

# Compute the median composite image using numpy
median_image = np.median(image_data, axis=0)

# Write the median composite image to a new geotif file
driver = gdal.GetDriverByName('GTiff')
output_dataset = driver.Create('path/to/output/image.tif', 
                                dataset.RasterXSize, 
                                dataset.RasterYSize, 
                                1, 
                                gdal.GDT_Float32)
output_dataset.SetGeoTransform(dataset.GetGeoTransform())
output_dataset.SetProjection(dataset.GetProjection())
output_dataset.GetRasterBand(1).WriteArray(median_image)
output_dataset.FlushCache()
```
