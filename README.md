# Data cube for Haskell University Workshop

<div style="width: 640px; height: 480px; margin: 10px; position: relative;"><iframe allowfullscreen frameborder="0" style="width:640px; height:480px" src="https://lucid.app/documents/embedded/7211b957-a462-4ea9-8894-64edf99dfef2" id="MGgkCvlPb5.R"></iframe></div> "title")

## 1. Create Area of Interest (Python code)

```python
# Import necessary modules first
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

### NLCD: National Land Cover Database 2019, Landcover & Imperviousness (NLCD2019)

Download from https://www.mrlc.gov/data?f%5B0%5D=region%3Aconus&f%5B1%5D=year%3A2019

### SRTM: 

## Create a median composite image

```python
import gdal
import numpy as np
import glob

# Get a list of all the geotif images in the specified directory
image_list = glob.glob('path/to/images/*.tif')

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