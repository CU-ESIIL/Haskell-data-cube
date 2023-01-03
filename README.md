# Data cube for Haskell University Workshop

![ESIIL_datacube](https://user-images.githubusercontent.com/67020853/210284955-b9cc208e-f736-4bc8-bec2-abc4fd2471b1.png)
Figure 1. Workflow for multi-source, satellite data processing

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
![Slide1](https://user-images.githubusercontent.com/67020853/210285112-724693da-193c-40b7-9e6d-c40d8423a8a4.GIF)
Figure 2. Area Of Interest (AOI) for Haskell University Workshop

## 2. Download data

### 2.1 NLCD: National Land Cover Database 2019, Landcover & Imperviousness (NLCD2019)

Download from https://www.mrlc.gov/data?f%5B0%5D=region%3Aconus&f%5B1%5D=year%3A2019

### 2.2 SRTM: NASA Shuttle Radar Topography Mission Global 1 arc second V003

Download from https://search.earthdata.nasa.gov/search

*Use the geojson file generated in section 1 for search

### 2.3 ECOSTRESS: ECOSTRESS Water Use Efficiency Daily L4 Global 70m V001 & ECOSTRESS Geolocation Daily L1B Global 70m V001

Download from https://search.earthdata.nasa.gov/search?q=ECOSTRESS

*Use the shape file generated in section 1 for search

**Restrict the L4 data search by year, season, and hours of the day

***Search and download L1B data by using the relative L4 data orbits 

### 2.4 GEDI: GEDI L1B Geolocated Waveform Data Global Footprint Level V002 & GEDI L2A Elevation and Height Metrics Data Global Footprint Level V002

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
## 3. Transform .img to .shp (NLCD) (Python code)
```python=
import gdal
import ogr

# Open the image file
image_file = gdal.Open('input.img')

# Get the image information and create the output file
driver = ogr.GetDriverByName('ESRI Shapefile')
image_info = image_file.GetGeoTransform()
output_file = driver.CreateDataSource('output.shp')
layer = output_file.CreateLayer('layer', image_info[1], ogr.wkbPoint)

# Loop through the image and create points for each pixel
for x in range(image_file.RasterXSize):
    for y in range(image_file.RasterYSize):
        # Get the pixel value
        pixel_value = image_file.GetRasterBand(1).ReadAsArray(x, y, 1, 1)[0][0]
        
        # Create a point for the pixel if it is not nodata
        if pixel_value != image_file.GetRasterBand(1).GetNoDataValue():
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(image_info[0] + x * image_info[1] + y * image_info[2], image_info[3] + x * image_info[4] + y * image_info[5])

            # Create a feature for the point and add it to the layer
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetGeometry(point)
            layer.CreateFeature(feature)

# Close the files
output_file.Destroy()
image_file = None
```
## 4. Reproject layer (NLCD) (Python code)
```python=
import os
from osgeo import ogr

# Set the input and output filenames
input_filename = 'input.shp'
output_filename = 'output.shp'

# Set the input and output spatial reference systems
input_srs = ogr.osr.SpatialReference()
input_srs.ImportFromEPSG(23700)  # Albert Conic projection
output_srs = ogr.osr.SpatialReference()
output_srs.ImportFromEPSG(4326)  # WGS84 projection

# Set the transformation
transform = ogr.osr.CoordinateTransformation(input_srs, output_srs)

# Open the input file
input_file = ogr.Open(input_filename)
input_layer = input_file.GetLayer()

# Create the output file
if os.path.exists(output_filename):
    os.remove(output_filename)
output_file = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(output_filename)
output_layer = output_file.CreateLayer(output_filename, output_srs, ogr.wkbPolygon)

# Add fields to the output file
for i in range(input_layer.GetLayerDefn().GetFieldCount()):
    field_defn = input_layer.GetLayerDefn().GetFieldDefn(i)
    output_layer.CreateField(field_defn)

# Add features to the output file
for feature in input_layer:
    # Transform the geometry
    geometry = feature.GetGeometryRef()
    geometry.Transform(transform)

    # Create a new feature
    output_feature = ogr.Feature(output_layer.GetLayerDefn())

    # Set the geometry and attributes
    output_feature.SetGeometry(geometry)
    for i in range(input_layer.GetLayerDefn().GetFieldCount()):
        output_feature.SetField(input_layer.GetLayerDefn().GetFieldDefn(i).GetNameRef(), feature.GetField(i))

    # Add the feature to the output file
    output_layer.CreateFeature(output_feature)

# Close the files
input_file = None
output_file = None
```

## 5. Mosaic scenes (SRTM) (Python code)

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
## 6. Transform .h5 to .tif (ECOSTRESS) (Python code)

Use script available at NASA EarthData Bitbucket: https://git.earthdata.nasa.gov/projects/LPDUR/repos/ecostress_swath2grid/browse

After setting-up the conda environment, run it in terminal as follows:

```python
> python ECOSTRESS_swath2grid.py --proj GEO --dir C:\Users\ECOSTRESS\
```
*If running the script in Windows/Git bash, type --dir path using forward slash

## 7. Create composite image (ECOSTRESS WUE) (Python code)

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

## 10. Mask data (Python code)

### 10.1 GeoTIFF data
```python=
import rasterio
from rasterio.mask import mask
from osgeo import ogr

# Open the TIF image
with rasterio.open("SRTM_Haskell.tif") as src:
    # Read the image data and metadata
    image_data, image_meta = src.read(masked=True)

# Open the shapefile
driver = ogr.GetDriverByName("ESRI Shapefile")
shapefile_ds = driver.Open("Haskell_wgs84.shp", 0)
shapefile_layer = shapefile_ds.GetLayer()

# Create an empty list to store the geometry objects
geoms = []

# Iterate through the features in the shapefile and extract the geometry
for feature in shapefile_layer:
    geom = feature.GetGeometryRef()
    geoms.append(geom)

# Use the rasterio mask function to clip the image to the shapefile boundary
clip_data, clip_meta = mask(src, geoms, crop=True)

# Save the clipped image
with rasterio.open("clipped_image.tif", "w", **clip_meta) as dst:
    dst.write(clip_data)

# Close the shapefile
shapefile_ds.Destroy()
```
### 10.2 Shapefile data
```python=
from osgeo import ogr

# Open the input shapefile
driver = ogr.GetDriverByName("ESRI Shapefile")
datasource = driver.Open("GEDI_merged.shp", 0)
layer = datasource.GetLayer()

# Open the boundary shapefile
boundary_ds = ogr.Open("Haskell_wgs84.shp")
boundary_layer = boundary_ds.GetLayer()
boundary_feature = boundary_layer.GetNextFeature()
boundary_geometry = boundary_feature.GetGeometryRef()

# Create a spatial filter based on the boundary geometry
layer.SetSpatialFilter(boundary_geometry)

# Write the filtered features to a new shapefile
out_ds = driver.CreateDataSource("clipped_shapefile.shp")
out_layer = out_ds.CreateLayer("clipped_shapefile", geom_type=layer.GetGeomType())
out_layer.CreateFields(layer.schema)
for feature in layer:
    out_layer.CreateFeature(feature)

# Close the shapefiles
datasource.Destroy()
boundary_ds.Destroy()
out_ds.Destroy()
```
![Slide2](https://user-images.githubusercontent.com/67020853/210285277-7b587bca-144c-49c7-addd-16c0a960acd7.GIF)
Figure 3. Land Cover (NLCD, Landsat-derived) from Haskell University Workshop AOI

![Slide3](https://user-images.githubusercontent.com/67020853/210285375-9566b3b9-92ea-4179-a14f-055eecea3f16.GIF)
Figure 4. Digital Elevation Model (SRTM) from Haskell University Workshop AOI

![Slide5](https://user-images.githubusercontent.com/67020853/210285593-c377495f-648e-44fd-a704-92af52f7cd7c.GIF)
Figure 5. Water Use Efficiency (ECOSTRESS L4) from Haskell University Workshop AOI 

![Slide6](https://user-images.githubusercontent.com/67020853/210285778-2ebd6436-cc92-49e5-9e97-ed2621173da8.GIF)
Figure 6. GEDI laser beams from Haskell University Workshop AOI 

![GEDI_samples_Haskell_slide10_500dpi](https://user-images.githubusercontent.com/67020853/210286109-3b1b7057-61a6-49ad-853a-a1668e2ba5bf.gif)
Figure 7. GEDI samples: pasture (yellow) and forest (cyan)

![WAVEFORMS](https://user-images.githubusercontent.com/67020853/210286536-d522f70a-7242-45df-886a-a85cf8111413.gif)
Figure 8. Pasture (A) and Forest (B) vertical profiles from Figure 7 GEDI samples. 