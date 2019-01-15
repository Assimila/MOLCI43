#!/bin/env python

from IPython import embed

import sys
import os
import numpy as np

HomeDir = os.path.expanduser('~')
SrcDir = os.path.join( HomeDir, 'Multiply/src/BRDF/QL')

from ScaleImage import ScaleImage

import osgeo.gdal as gdal
from osgeo.gdalconst import *

PriorFile = sys.argv[1]
tile = sys.argv[2]
mask_fname = sys.argv[3]
d_mask = gdal.Open(mask_fname)
mask = d_mask.ReadAsArray()
mask = np.where((mask >= 1) & (mask <=3), True, False)

# False color composite, RGB 721
# R: Band 19 Mean - band: 7 Parameter f0
# G: Band 4  Mean - band: 2 Parameter f0
# B: Band 1  Mean - band: 1 Parameter f0

dataset = gdal.Open( PriorFile )

Projection = dataset.GetProjection()
GeoTransform = dataset.GetGeoTransform()

f0_b7 = dataset.GetRasterBand( 13 ).ReadAsArray() * mask
f0_b7_scaled = ScaleImage( f0_b7 )

f0_b2 = dataset.GetRasterBand( 4 ).ReadAsArray() * mask
f0_b2_scaled = ScaleImage( f0_b2 )

f0_b1 = dataset.GetRasterBand( 1 ).ReadAsArray() * mask
f0_b1_scaled = ScaleImage( f0_b1 )

# Save scaled image to a GeoTiff file
format = "GTiff"
driver = gdal.GetDriverByName(format)

NumberOfBands = 1
#cols, rows = f0_b1_scaled.shape[1], f0_b1_scaled.shape[0]
cols, rows = f0_b1.shape[1], f0_b1.shape[0]

new_dataset = driver.Create( 'Scaled_b1.' + tile + '.tif', cols, rows, NumberOfBands, GDT_Byte )
new_dataset.GetRasterBand(1).WriteArray(f0_b1_scaled)
new_dataset.SetProjection( Projection )
new_dataset.SetGeoTransform( GeoTransform )
new_dataset = None

new_dataset = driver.Create( 'Scaled_b2.' + tile + '.tif', cols, rows, NumberOfBands, GDT_Byte )
new_dataset.GetRasterBand(1).WriteArray(f0_b2_scaled)
new_dataset.SetProjection( Projection )
new_dataset.SetGeoTransform( GeoTransform )
new_dataset = None

new_dataset = driver.Create( 'Scaled_b7.' + tile + '.tif', cols, rows, NumberOfBands, GDT_Byte )
new_dataset.GetRasterBand(1).WriteArray(f0_b7_scaled)
new_dataset.SetProjection( Projection )
new_dataset.SetGeoTransform( GeoTransform )
new_dataset = None

# True color composie, RGB 143
# R: Band 1  Mean - band: 1 Parameter f0
# G: Band 10 Mean - band: 4 Parameter f0
# G: Band 7  Mean - band: 3 Parameter f0

f0_b4 = dataset.GetRasterBand( 10 ).ReadAsArray() * mask
f0_b4_scaled = ScaleImage( f0_b4 )

f0_b3 = dataset.GetRasterBand( 7 ).ReadAsArray() * mask
f0_b3_scaled = ScaleImage( f0_b3 )

# Save scaled image to a GeoTiff file
format = "GTiff"
driver = gdal.GetDriverByName(format)

new_dataset = driver.Create( 'Scaled_b4.' + tile + '.tif', cols, rows, NumberOfBands, GDT_Byte )
new_dataset.SetProjection( Projection )
new_dataset.SetGeoTransform( GeoTransform )
new_dataset.GetRasterBand(1).WriteArray(f0_b4_scaled)
new_dataset = None

new_dataset = driver.Create( 'Scaled_b3.' + tile + '.tif', cols, rows, NumberOfBands, GDT_Byte )
new_dataset.SetProjection( Projection )
new_dataset.SetGeoTransform( GeoTransform )
new_dataset.GetRasterBand(1).WriteArray(f0_b3_scaled)

new_dataset = None

