#!/bin/env python

from IPython import embed

import sys
import os

HomeDir = os.path.expanduser('~')
SrcDir = os.path.join( HomeDir, 'Multiply/src/BRDF/QL')

from ScaleImage import ScaleImage

import osgeo.gdal as gdal
from osgeo.gdalconst import *

fname = sys.argv[1]
band = int(sys.argv[2])
dataset = gdal.Open(fname)

Projection = dataset.GetProjection()
GeoTransform = dataset.GetGeoTransform()

b = dataset.GetRasterBand(1).ReadAsArray()
b_scaled = ScaleImage(b)

# Save scaled image to a GeoTiff file
format = "GTiff"
driver = gdal.GetDriverByName(format)

NumberOfBands = 1
rows, cols = b_scaled.shape

new_dataset = driver.Create( 'Scaled_b%d.tif' % band, cols, rows, NumberOfBands, GDT_Byte )
new_dataset.GetRasterBand(1).WriteArray(b_scaled)
new_dataset.SetProjection( Projection )
new_dataset.SetGeoTransform( GeoTransform )
new_dataset = None
