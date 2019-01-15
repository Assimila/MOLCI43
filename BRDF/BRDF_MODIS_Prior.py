#!/bin/env python

import glob
import os
import sys

import numpy as np
import osgeo.gdal as gdal
from osgeo.gdalconst import *

def GetProjectionParams( File ):
    """
    Purpose
    Get projection parameters from the MCD43A1 product
    """
    dataset = gdal.Open( File )
    SubDatasets = dataset.GetSubDatasets()

    SubDataset = gdal.Open( SubDatasets[0][0] )
    Projection = SubDataset.GetProjection()
    GeoTransform = SubDataset.GetGeoTransform()

    return Projection, GeoTransform

def IncrementSamples(AccWeightedParams, AccWeight, Parameters, Uncertainties):
    """
    Author
        Gerardo Lopez Saldana, G.LopezSaldana@reading.ac.uk

    Purpose
    Accumulate information applying a weighted samples factor 

    Inputs
    AccWeightedParams: Float[rows, cols, NumberOfBands, NumberOfParameters]: Accumulator array to store the weigthed contribution of Params
    AccWeight: Float[rows, cols, NumberOfBands]: Array of the sum of weighting factors
    Parameters: Float[rows, cols, NumberOfBands, NumberOfParameters]: Array containg the BRDF parameters (iso,vol,geo)
    Uncertainties: Float[rows, cols, NumberOfBands]: Array containing the uncertainty information deraived from the QA flag

    Outputs
    Acumulator
    """
    rows, cols, NumberOfBands, NumberOfParameters = Parameters.shape

    for i in range(NumberOfBands):
        # Create the array for the weigthing factors
        Weight = np.where(Uncertainties[:,:,i] > 0, Uncertainties[:,:,i], 0.0)

        for j in range(NumberOfParameters):
            AccWeightedParams[:,:,i,j] += Parameters[:,:,i,j] * Weight

        AccWeight[:,:,i] += Weight

    return AccWeightedParams, AccWeight

def GetParameters(ParamsFile, QualityFile, Bands, NumberOfParameters, RelativeUncert, ScaleFactor, ProcessSnow = 0):
    """
    Author
        Gerardo Lopez Saldana, G.LopezSaldana@reading.ac.uk

    Purpose
        Method to extract BRDF parameters from MCD43A1  SDS

    Inputs
        ParamsFile: String: MCD4A2 File name (absolute path)
        QualityFile: String: MCD4A1 File name (absolute path)
        Bands: Long[]: Array containing the bands to extract (1,2 == Red, NIR)
        NumberOfParameters: Long[]: Array containing the number of BRDF parameters
        RelativeUncert: Float[]: 5 elements array containing the relative uncertainty for each QA value

    Outputs
        Parameters: Array containg the BRDF parameters (iso,vol,geo) for Bands
        Uncertainties: Array containing the uncertainty information deraived from the QA flag
    """

    FillValue = 32767
    NumberOfBands = Bands.shape[0]

    # Get dimensions
    rows, cols = GetDimSubDataset( ParamsFile )

    Parameters = np.zeros((rows, cols, NumberOfBands, NumberOfParameters), np.float32)
    Uncertainties = np.zeros((rows, cols, NumberOfBands), np.float32)

    # Get Snow
    # 1     Snow albedo retrieved
    # 0     Snow-free albedo retrieved
    # 255   Fill Value
    print "Reading Snow QA:", QualityFile
    SubDatasetName = 'HDF4_EOS:EOS_GRID:"' + QualityFile + '":MOD_Grid_BRDF:Snow_BRDF_Albedo'
    SubDataset = gdal.Open(SubDatasetName, GA_ReadOnly)
    SnowQA = SubDataset.GetRasterBand(1).ReadAsArray()
    if ProcessSnow == 0:
        SnowQA = np.where( SnowQA == 0, 1, 0)
    else:
        SnowQA = np.where( SnowQA == 1, 1, 0)

    # Load BRDF parameters
    print "Reading BRDF parameters..."
    for Band in range( Bands.shape[0] ):
        SubDatasetName = 'HDF4_EOS:EOS_GRID:"' + ParamsFile + '":MOD_Grid_BRDF:BRDF_Albedo_Parameters_Band' + str( Bands[Band] )
        print SubDatasetName 
        SubDataset = gdal.Open(SubDatasetName, GA_ReadOnly)

        for Parameter in range(NumberOfParameters):
            print "Getting BRDF parameter", Parameter
            Parameters[:,:,Band,Parameter] = SubDataset.GetRasterBand( Parameter + 1 ).ReadAsArray()

            # Snow mask
            Parameters[:,:,Band,Parameter] = Parameters[:,:,Band,Parameter] * SnowQA

            # Filter out fill values
            Parameters[:,:,Band,Parameter] = np.where(Parameters[:,:,Band,Parameter] == FillValue, 0.,
                                                      Parameters[:,:,Band,Parameter] * ScaleFactor )

    # Get QA
    print "Reading QA:", QualityFile
    for Band in range( Bands.shape[0] ):
        SubDatasetName = 'HDF4_EOS:EOS_GRID:"' + QualityFile + '":MOD_Grid_BRDF:BRDF_Albedo_Band_Quality_Band' + str( Bands[Band] )
        SubDataset = gdal.Open(SubDatasetName, GA_ReadOnly)
        QA = SubDataset.GetRasterBand(1).ReadAsArray()

        # https://ladsweb.nascom.nasa.gov/api/v1/filespec/collection=6&product=MCD43A2
        # BRDF_Albedo_Band_Quality_BandN ( N is 1 to 7 )> 
        # 0 = best quality, full inversion (WoDs, RMSE majority good)
        # 1 = good quality, full inversion (also including the cases that no clear sky
        #     observations over the day of interest or the Solar Zenith Angle is too 
        #     large even WoDs, RMSE majority good)
        # 2 = Magnitude inversion (numobs >=7)
        # 3 = Magnitude inversion (numobs >=2&<7)
        # 4 = Fill value

        QA_flags = np.array( [ 0,1,2,3 ] )

        for i, QA_flag in enumerate( QA_flags ) :
            indices = np.where( QA == QA_flag )
            Uncertainties[ indices[0], indices[1], Band ] = RelativeUncert[ i ]

        Uncertainties[:,:,Band] =  Uncertainties[:,:,Band] * SnowQA        

    SubDataset = None
    return Parameters, Uncertainties

def GetFileList(DataDir, Product, DoY, Collection):
    # DataDir

    FileList = glob.glob( os.path.join( DataDir,
        Product, "%s.A20??%s.h??v??.%s.*.hdf" % ( Product, DoY, Collection ) ) )
    FileList.sort()

    Year = np.zeros((len(FileList)), np.int16)

    i = 0
    for File in FileList:
        # Get Year from filename
        YearOfObservation = os.path.basename(File).split('.')[1][1:5]
        Year[i] = YearOfObservation

        i += 1

    return FileList, Year


def GetDimSubDataset(File):
    # Open first subdataset from MCD43A?
    Dataset = gdal.Open(File, GA_ReadOnly)
    SubDataset = Dataset.GetSubDatasets()[0][0]

    dataset = gdal.Open(SubDataset, GA_ReadOnly)
    rows, cols = dataset.RasterYSize, dataset.RasterXSize
    dataset = None

    return rows, cols

#--------------------------------------------------------------------------------#
from IPython import embed

#HomeDir = os.path.expanduser('~')
#DataDir = os.path.join( HomeDir , 'MELODIES/data/MODIS' )
DataDir = "/data/MODIS/h09v07"

DoY = sys.argv[1]
DoY = "%03d" % ( int( DoY ) )

ProcessSnow = int( sys.argv[2] )
Collection = '006'

# BRDF parameters
Product = 'MCD43A1'
FileListA1, Year = GetFileList( DataDir, Product, DoY, Collection )

# BRDF Quality
Product = 'MCD43A2'
FileListA2, Year = GetFileList( DataDir, Product, DoY, Collection )

if len(FileListA1) <> len(FileListA2):
    print 'File lists are inconsistent.'
    exit(-1)

# From the first file get dimensions and projection information
rows, cols = GetDimSubDataset( FileListA1[0] )
Projection, GeoTransform =  GetProjectionParams( FileListA1[0] )

# Array of bands to process
Bands = np.array( [1,2,3,4,5,6,7] )
NumberOfBands = len( Bands )
# BRDF parameters, iso, vol, geo
NumberOfParameters = 3

ScaleFactor = 0.001

# Relative uncertainty of 4 quality values
# https://ladsweb.nascom.nasa.gov/api/v1/filespec/collection=6&product=MCD43A2
BRDF_Albedo_Quality = np.arange(4)

# http://en.wikipedia.org/wiki/Golden_ratio
GoldenMean = 0.618034
RelativeUncert = GoldenMean ** BRDF_Albedo_Quality

WeightedMean = np.zeros((rows, cols, NumberOfBands, NumberOfParameters), np.float32)
WeightedVariance = np.zeros((rows, cols, NumberOfBands, NumberOfParameters), np.float32)
AccWeightedParams = np.zeros((rows, cols, NumberOfBands, NumberOfParameters), np.float32)
AccWeight = np.zeros((rows, cols,NumberOfBands), np.float32)

print "Computing the weigthed mean..."
for A1file, A2file in zip( FileListA1, FileListA2 ):
    Parameters, Uncertainties = GetParameters(A1file, A2file, Bands, NumberOfParameters,
                                              RelativeUncert, ScaleFactor, ProcessSnow)
    AccWeightedParams, AccWeight = IncrementSamples(AccWeightedParams, AccWeight, Parameters, Uncertainties)

# Compute the weighted mean
for i in range(NumberOfBands):
    for j in range(NumberOfParameters):
        WeightedMean[:,:,i,j] = np.where(AccWeight[:,:,i] > 0.0, AccWeightedParams[:,:,i,j] / AccWeight[:,:,i], 0.0)

print "Computing the weigthed variance..."
# Compute the weigthed variance
for A1file, A2file in zip( FileListA1, FileListA2 ):
    Parameters, Uncertainties = GetParameters(A1file, A2file, Bands, NumberOfParameters,
                                              RelativeUncert, ScaleFactor, ProcessSnow)
    for i in range(NumberOfBands):
        for j in range(NumberOfParameters):
            tmpWeightedVariance = Uncertainties[:,:,i] * np.power(Parameters[:,:,i,j] - WeightedMean[:,:,i,j], 2)
            WeightedVariance[:,:,i,j] += tmpWeightedVariance

for i in range(NumberOfBands):
    for j in range(NumberOfParameters):
        WeightedVariance[:,:,i,j] = np.where(AccWeight[:,:,i] > 0.0, WeightedVariance[:,:,i,j] / AccWeight[:,:,i], 0.0)

print "Writing results to a file..."
format = "GTiff"
driver = gdal.GetDriverByName(format)
if ProcessSnow == 0:
    OutputFilename = 'MCD43A1.Prior.' + DoY +  '.img'
else:
    OutputFilename = 'MCD43A1.SnowPrior.' + DoY +  '.img'

new_dataset = driver.Create( OutputFilename, cols, rows, \
                             ((NumberOfBands*NumberOfParameters)*2)+NumberOfBands, 
                             GDT_Float32, ['COMPRESS=DEFLATE', 'BIGTIFF=YES',
                         'PREDICTOR=1', 'TILED=YES'] )

k = 1
for i in range(NumberOfBands):
    for j in range(NumberOfParameters):
        new_dataset.GetRasterBand(k).WriteArray(WeightedMean[:,:,i,j])
        new_dataset.GetRasterBand(k).SetDescription("Mean - band: " + str(Bands[i]) + " Parameter f" + str(j))

        new_dataset.GetRasterBand((NumberOfBands*NumberOfParameters)+k).SetDescription("VAR - band: " + str(Bands[i]) + " Parameter f" + str(j))
        new_dataset.GetRasterBand((NumberOfBands*NumberOfParameters)+k).WriteArray(WeightedVariance[:,:,i,j])
        k += 1

    # Write the weigthed number of samples
    new_dataset.GetRasterBand(((NumberOfBands*NumberOfParameters)*2) + i + 1).SetDescription("Weighted Number Of Samples Band: " + str(Bands[i]))
    new_dataset.GetRasterBand(((NumberOfBands*NumberOfParameters)*2) + i + 1).WriteArray(AccWeight[:,:,i])

new_dataset.SetProjection( Projection )
new_dataset.SetGeoTransform( GeoTransform )

new_dataset = None


