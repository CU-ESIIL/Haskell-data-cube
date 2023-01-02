# Data cube for Haskell University Workshop

![ESIIL_datacube](https://user-images.githubusercontent.com/67020853/210265589-598490f7-113b-4c9a-8cbf-48b67ab0a5c5.png)

## 1. Create Area of Interest (Python code)

```python
# Import 
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

### SRTM: NASA Shuttle Radar Topography Mission Global 1 arc second V003

Download from https://search.earthdata.nasa.gov/search

*Use the geojson file generated in section 1 for search

### ECOSTRESS: ECOSTRESS Water Use Efficiency Daily L4 Global 70m V001 & ECOSTRESS Geolocation Daily L1B Global 70m V001

Download from https://search.earthdata.nasa.gov/search?q=ECOSTRESS

*Use the shape file generated in section 1 for search

**Restrict the L4 data search by year, season, and hours of the day

***Search and download L1B data by using the relative L4 data orbits 

### GEDI: GEDI L1B Geolocated Waveform Data Global Footprint Level V002 & GEDI L2A Elevation and Height Metrics Data Global Footprint Level V002

Download using the rGEDI package: https://github.com/carlos-alberto-silva/rGEDI

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