# -*- coding: utf-8 -*-
"""
Purpose
  The MOD09GA script performs:
    - Cloud screening based on the MOD09GA QA flags
    - Sets an absolute uncertainty based on the AOD quantities
    - Computes the Ross-Thick and Li-Sparse kernels

  To process MOD09GA for 2015:
  python MOD09GA.py MOD09GA 2015

Inputs
  arg 1 - Product e.g. MOD09GA or MYD09GA
  arg 2 - Year e.g. 2015

  The MOD09GA scripts requites a configuration file so it can
  know the location of the input data and format, the output
  directory and format, tile, etc. Default is .BEIS_LC.cfg 

Outputs
  For every input file two files will be produced, e.g.:
    - MOD09GA.2014347.h17v03.tif
      Containing all first seven MODIS reflectance bands with
      screened of clouds and associated uncertainties for year
      2014, julian day 347, tile h17v03 in GeoTiff format.
    - MOD09GA.2014347.h17v03.kernels.tif
      The Ross-Thick and Li-Sparse kernels for the particular
      viewing and illumination conditions of the scene
"""

__author__ = "Gerardo López Saldaña"
__version__ = "0.2 (19.04.2017)"
__email__ = "gerardo.lopezsaldana@assimila.eu"

import sys
import os
import glob
import configparser as ConfigParser

try:
    import numpy as np
except ImportError:
    print('Numpy is not installed.')
    exit(-1)

try:
    import osgeo.gdal as gdal
    from osgeo.gdalconst import *
    gdal.UseExceptions()
except ImportError:
    print('GDAL is not installed.')
    exit(-1)

from IPython import embed

if len ( sys.argv ) != 3:
    print(__doc__)
    exit( -1 )

# Set directory structure
# Get configuration
config = ConfigParser.ConfigParser()
HomeDir = os.path.expanduser('~')
config.read( os.path.join( HomeDir, 'Multiply/src/BRDF/.Multiply_BRDF.cfg' ) )

Project = config.get( 'Project', 'Project' )
cDataDir = config.get( 'Directory', 'DataDir' )
cSrcDir = config.get( 'Directory', 'SrcDir' )
cScreenedDataDir = config.get( 'Directory', 'ScreenedDataDir' )
Sensor = config.get( 'Data', 'Sensor' )
#Tile = config.get( 'Data', 'Tile' )
Tile = 'h17v05'
#InputFormat = config.get( 'Data', 'InputFormat' )
InputFormat = 'hdf'
#NumberOfReflBands = config.getint( 'Data', 'NumberOfReflBands' )
NumberOfReflBands = 7

Product = sys.argv[1]
Year = sys.argv[2]

# $HOME/PROJECT/data/MODIS/TILE/PRODUCT/YEAR
DataDir = os.path.join( HomeDir , Project, cDataDir, Sensor, Tile, Product, Year )
OutputDir = os.path.join( DataDir, cScreenedDataDir )

SrcDir = os.path.join( HomeDir , Project, cSrcDir, Sensor )
sys.path.append( SrcDir )

Files = glob.glob( os.path.join( DataDir, '*' + InputFormat ) )
print(os.path.join( DataDir, '*' + InputFormat ))
Files.sort()

for File in Files:
    print(File)
    # Internal HDF file structure
    # HDF4_EOS:EOS_GRID:"FILE":MODIS_Grid_500m_2D:sur_refl_b0N_1
    for Band in range( NumberOfReflBands ):
        dataset = gdal.Open( 'HDF4_EOS:EOS_GRID:"%s":MODIS_Grid_500m_2D:sur_refl_b0%s_1' % ( File, Band + 1 ) )
        scale_factor = float( dataset.GetMetadataItem( 'scale_factor' ) )
        # Get projection information
        Projection = dataset.GetProjection()
        GeoTransform = dataset.GetGeoTransform()

        # If laterstack exists
        if 'layerstack' in locals():
            layerstack[ Band ] = dataset.ReadAsArray()
        else:
            rows, cols = dataset.RasterYSize, dataset.RasterXSize
            # Layerstack for reflectances
            layerstack = np.zeros( ( NumberOfReflBands, rows, cols ), np.int16 )
            layerstack[ Band ] = dataset.ReadAsArray()

        dataset = None

    # QA information
    # Bits are listed from the MSB (bit 31) to the LSB (bit 0):
    # Bit    Description
    # 2      cloud shadow;
    #        1 -- yes
    #        0 -- no
    # 0-1    cloud state;
    #        00 -- clear
    #        01 -- cloudy
    #        10 -- mixed
    #        11 -- not set, assumed clear

    # Internal HDF file structure
    # HDF4_EOS:EOS_GRID:"FILE":MODIS_Grid_1km_2D:state_1km_1
    dataset = gdal.Open( 'HDF4_EOS:EOS_GRID:"%s":MODIS_Grid_1km_2D:state_1km_1' % ( File ) )
    QA = dataset.ReadAsArray()
    # Downscale the 1km QA mask to 500m, hence a factor of 2
    ResamplingScaleFactor = 2
    ResamplingArray = np.ones( ( ResamplingScaleFactor, ResamplingScaleFactor ), np.int8 )

    QA = np.kron( QA , ResamplingArray )

    # Cloud state 00
    bit0, bit1 = 0, 1 
    Cloud = np.where( ( ( ( QA / np.power( 2, bit0 ) ) % 2 == 0 ) & 
                        ( ( QA / np.power( 2, bit1 ) ) % 2 == 0 )  ), 1 , 0 )
 
    # Bit 2 cloud shadow: 0 - no cloud shadow
    bit2 = 2
    CloudShadow = np.where( ( QA / np.power( 2, bit2 ) ) % 2 == 0, 1, 0)

    # Bits 3-5 land-water flag
    # 000 - shallow ocean
    # 110 - continental/moderate ocean
    # 111 - deep ocean
    bit3 = 3
    bit4 = 4
    bit5 = 5
    # 000 - shallow ocean
    LandWater =  np.where( ( ( QA / np.power( 2, bit3 ) ) % 2 == 0) & \
                           ( ( QA / np.power( 2, bit4 ) ) % 2 == 0) & \
                           ( ( QA / np.power( 2, bit5 ) ) % 2 == 0) , 0, 1)
    # 110 - continental/moderate ocean
    LandWater =  np.where( ( ( QA / np.power( 2, bit3 ) ) % 2 == 0) & \
                           ( ( QA / np.power( 2, bit4 ) ) % 2 == 1) & \
                           ( ( QA / np.power( 2, bit5 ) ) % 2 == 1) , 0, LandWater)
    # 111 - deep ocean
    LandWater =  np.where( ( ( QA / np.power( 2, bit3 ) ) % 2 == 1) & \
                           ( ( QA / np.power( 2, bit4 ) ) % 2 == 1) & \
                           ( ( QA / np.power( 2, bit5 ) ) % 2 == 1) , 0, LandWater)

    # Bit 10  internal cloud flag: 0 - no cloud
    bit10 = 10
    #InternalCloudFlag = np.where( ( QA / np.power( 2, bit10 ) ) % 2 == 0, 1, 0 )
    InternalCloudFlag = np.where( ( QA / np.power( 2, bit10 ) ) % 2 == 0, 0, 1)

    # Bit 12 MOD35 snow/ice flag: 1 - snow
    bit12 = 12
    #SnowIceFlag = np.where( ( QA / np.power( 2, bit12 ) ) % 2 == 1, 1, 0 )
    SnowIceFlag = np.where( ( QA / np.power( 2, bit12 ) ) % 2 == 1, 0, 1 )

    # Bit 13 Pixel is adjacent to cloud : 0 - no
    bit13 = 13
    #PixelAdjacentToCloud = np.where( ( QA / np.power( 2, bit13 ) ) % 2 == 0, 1, 0 )
    PixelAdjacentToCloud = np.where( ( QA / np.power( 2, bit13 ) ) % 2 == 0, 0, 1 )

    # Set uncert based aerosol QA
    # Bit 6-7 Aerosol quantity:
    # 00 - climatology
    # 01 - low
    # 10 - average
    # 11 - high
    Aerosols = { 1:0.01, 2:0.02, 3:0.03, 4:0.04 }
    bit6, bit7 = 6, 7

    UncerntAOD = np.array(4, np.float32)
    AerosolsQA = np.where( ((QA / np.power(2,bit6)) % 2 == 0) & ((QA / np.power(2,bit7)) % 2 == 0), Aerosols[1], 0.0 )
    AerosolsQA = np.where( ((QA / np.power(2,bit6)) % 2 == 0) & ((QA / np.power(2,bit7)) % 2 == 1), Aerosols[2], AerosolsQA)
    AerosolsQA = np.where( ((QA / np.power(2,bit6)) % 2 == 1) & ((QA / np.power(2,bit7)) % 2 == 0), Aerosols[3], AerosolsQA)
    AerosolsQA = np.where( ((QA / np.power(2,bit6)) % 2 == 1) & ((QA / np.power(2,bit7)) % 2 == 1), Aerosols[4], AerosolsQA)

    refl_SD_noise_estimates = {1:0.004, 2:0.015, 3:0.003, 4:0.004, 5:0.013, 6:0.010, 7:0.006 }

    # Save masked refl to a GeoTiff file
    format = "GTiff"
    driver = gdal.GetDriverByName(format)

    # fname e.g. MOD09GA.A2015235.h17v03.006.2015305215014.hdf
    fname = os.path.basename( File )
    Date = fname.split('.')[1][1::]
    fname = '%s.%s.%s.tif' % ( Product, Date, Tile ) 
    fname = os.path.join( OutputDir , fname )

    # Seven spectral bands plus uncert and snow mask
    new_dataset = driver.Create( fname, cols, rows, ( NumberOfReflBands * 2) + 1, 
                  GDT_Int16 , options=["COMPRESS=LZW", "INTERLEAVE=BAND", "TILED=YES"] )

    # Layerstack for reflectances uncert
    Uncert = np.zeros( ( NumberOfReflBands, rows, cols ), np.int16 ) 

    # Screening using QA
    for i in range( NumberOfReflBands ):
        layerstack[i] = layerstack[i] * Cloud * CloudShadow * LandWater * InternalCloudFlag * PixelAdjacentToCloud
        layerstack[i] = np.where( (layerstack[i] < 0) | ( layerstack[i] > scale_factor ), 0, layerstack[i] )
        new_dataset.GetRasterBand(i + 1).WriteArray( layerstack[i] )

        Uncert[i] = ( ( ( layerstack[i] / scale_factor ) * AerosolsQA ) + \
                    ( refl_SD_noise_estimates[i + 1] ) ) * scale_factor

        Uncert[i] = Uncert[i]  * Cloud * CloudShadow * LandWater * InternalCloudFlag * PixelAdjacentToCloud
        Uncert[i] = np.where( ( Uncert[i]  < 0 ) | ( Uncert[i] > scale_factor ), 0, Uncert[i] )
        new_dataset.GetRasterBand( i + 1 + NumberOfReflBands).WriteArray( Uncert[i] )

    SnowMask = SnowIceFlag * Cloud * CloudShadow * LandWater * InternalCloudFlag * PixelAdjacentToCloud
    new_dataset.GetRasterBand( (NumberOfReflBands * 2) + 1 ).WriteArray( SnowMask )

    new_dataset.SetProjection( Projection )
    new_dataset.SetGeoTransform( GeoTransform )

    new_dataset = None

    # Let's tyde up a bit before calculating the kernels
    del layerstack , Uncert

    #==================
    # Calculate Kernels
    #==================
    from kernels import Kernels

    # View Zenith Angles - VZA
    dataset = gdal.Open( 'HDF4_EOS:EOS_GRID:"%s":MODIS_Grid_1km_2D:SensorZenith_1' % ( File ) )
    scale_factor = float( dataset.GetMetadataItem( 'scale_factor' ) )
    VZA = dataset.ReadAsArray() * scale_factor 

    # Solar Zenith Angles - SZA
    dataset = gdal.Open( 'HDF4_EOS:EOS_GRID:"%s":MODIS_Grid_1km_2D:SolarZenith_1' % ( File ) )
    SZA = dataset.ReadAsArray() * scale_factor 

    # Sensor Azimuth Angles
    dataset = gdal.Open( 'HDF4_EOS:EOS_GRID:"%s":MODIS_Grid_1km_2D:SensorAzimuth_1' % ( File ) )
    VAA = dataset.ReadAsArray() * scale_factor 

    # Solar Zenith Angles
    dataset = gdal.Open( 'HDF4_EOS:EOS_GRID:"%s":MODIS_Grid_1km_2D:SolarZenith_1' % ( File ) )
    SAA = dataset.ReadAsArray() * scale_factor 

    # Calculate Relative Azimuth Angle as
    # Solar azimuth minus viewing azimuth
    RAA = SAA - VAA

    # Kernels receives vza, sza, raa
    print('Calculating Ross-Thick and Li-Sparse kernels...')
    kk = Kernels(VZA, SZA, RAA, \
             RossHS = False, RecipFlag = True, MODISSPARSE = True, \
             normalise = 1, doIntegrals = False, LiType = 'Sparse', RossType = 'Thick' )

    Ross_Thick = np.zeros( ( rows, cols ), np.float32 )
    Li_Sparse = np.zeros( ( rows, cols ), np.float32 )

    # Downscale kernels to match the reflectance pixel size
    ResamplingScaleFactor = 2
    ResamplingArray = np.ones( ( ResamplingScaleFactor, ResamplingScaleFactor ), np.float32 )

    Ross_Thick = kk.Ross.reshape ( ( int(rows / ResamplingScaleFactor) , int(cols / ResamplingScaleFactor) ) )
    Ross_Thick = np.kron( Ross_Thick , ResamplingArray )
    Ross_Thick = Ross_Thick * Cloud * CloudShadow * LandWater * InternalCloudFlag * PixelAdjacentToCloud

    Li_Sparse = kk.Li.reshape( ( int(rows / ResamplingScaleFactor) , int(cols / ResamplingScaleFactor) ) )
    Li_Sparse = np.kron( Li_Sparse , ResamplingArray )
    Li_Sparse = Li_Sparse * Cloud * CloudShadow * LandWater * InternalCloudFlag * PixelAdjacentToCloud

    kk = None

    # Save kernels in a different file to keep all 32bit floating point information
    fname = '%s.%s.%s.kernels.tif' % ( Product, Date, Tile )
    fname = os.path.join( OutputDir , fname )

    new_dataset = driver.Create( fname, cols, rows, 2, GDT_Float32,
                      options = ["COMPRESS=LZW", "INTERLEAVE=BAND", "TILED=YES"] )

    new_dataset.GetRasterBand(1).WriteArray( Ross_Thick )
    new_dataset.GetRasterBand(2).WriteArray( Li_Sparse )

    new_dataset.SetProjection( Projection )
    new_dataset.SetGeoTransform( GeoTransform )

    new_dataset = None

