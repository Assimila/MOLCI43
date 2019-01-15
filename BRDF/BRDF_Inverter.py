# -*- coding: utf-8 -*-

"""
Purpose
  The BRDF_Inverter.py script:
    - Derives the BRDF model parameters for a particular
      day of the yeaar (DoY) using all available daily surface
      reflectance observations within a 32-day window
      to invert the BRDF model. The DoY is centered at the
      32-day window. The BRDF model parameters are used to
      compute reflectance under any valid viewing and
      illumination conditions.

  To derive the BRDF model parameters for January 1, 2011:
  python BRDF_Inverter.py 1 2011

Inputs
  arg 1 - The julian day to process, e.g. 1
  arg 2 - Year e.g. 2011

  The BRDF_Inverter.py scripts requites a configuration file
  so it can know the location of the input data and format, 
  the output directory and format, tile, etc. 
  Default is .BEIS_LC.cfg 

Outputs
    - A compressed GeoTiff file: BRDF_Parameters.YYYYddd.tif
    where YYYY is the year selected and ddd the julian day
    using three digits.
    The files has 44 bands.
      - First 21 are the BRDF model parameters for each of the
      MODIS bands 1 to 7.
      - Bands 22 to 42 the model parameters uncertainties.
      - Band 43 the number of samples used.
      - Band 44 the goodness of fit.
"""

__author__ = "Gerardo López Saldaña"
__version__ = "0.2 (20.04.2017)"
__email__ = "gerardo.lopezsaldana@assimila.eu"

import glob
import os
import sys
import time
import configparser as ConfigParser

from multiprocessing import Process, Array
import numpy as np
import datetime as DT
from matplotlib.dates import date2num
import math
from scipy.linalg import lstsq, inv
from scipy.optimize import nnls

try:
    import sharedmem as shm
except ImportError:
    print('Numpy/SharedMemory is not installed.')

try:
    import statsmodels.api as sm
except ImportError:
    print('StatsModel lib is not installed')
    exit(-1)

try:
    import osgeo.gdal as gdal
    from osgeo.gdalconst import *
    gdal.UseExceptions()
except ImportError:
    print('GDAL is not installed.')
    exit(-1)

def NormalDistribution(mu,sigma):
    def f(x):
        z = 1.0 * (x - mu) / sigma
        e = math.e ** (-0.5 *z ** 2)
        C = math.sqrt(2 * math.pi) * sigma
        return 1.0 * e / C
    return f

def AssessObservations( NIR_profile, K_profiles ):
    indices = np.where( NIR_profile > 0.0 )
    if indices[0].shape[0] < 5:
        valid_obs = True
        return valid_obs

    K = np.ones( [ 3, indices[0].shape[0] ] )
    K[1,:] = K_profiles[ indices, 1 ][0]
    K[2,:] = K_profiles[ indices, 2 ][0]

    K = np.matrix( K )
    M_ = K * K.T
    MI = M_.I

    R = NIR_profile[ [indices] ]
    V_ = K * R.T
    # Get parameters
    P_ = ( MI * V_ ).T

    # Get RMSE
    FWD = P_ * K
    d = FWD - R
    e = np.dot( d[0], d[0].T )
    rmse = np.sqrt( e / indices[0].shape[0] )

    # Get BHR
    BHR = P_[0,0] + ( P_[0,1] * 0.189184 ) + ( P_[0,2]* -1.377622 )
    if ( BHR > ( R.min() - rmse ) ) and ( BHR < ( R.max() + rmse ) ):
        valid_obs = True
    else:
        valid_obs = False

    return valid_obs

def GetOutlierIndices( NIR_profile, K_profiles, error_threshold = 0.3 ):

    indices = np.where( NIR_profile > 0.0 )
    if indices[0].shape[0] < 5:
        outlier_indices = np.array([])
        return outlier_indices

    K = np.ones( [ 3, indices[0].shape[0] ] )
    K[1,:] = K_profiles[ indices, 1 ][0]
    K[2,:] = K_profiles[ indices, 2 ][0]

    K = np.matrix( K )
    M_ = K * K.T
    MI = M_.I

    R = NIR_profile[ [indices] ]
    V_ = K * R.T
    # Get parameters
    P_ = ( MI * V_ ).T

    # Get relative error from forwarded model
    FWD = P_ * K
    d = FWD - R
    rel_err = np.abs( d / R )

    # Get RMSE
    e = np.dot( d[0], d[0].T )
    rmse = np.sqrt( e / indices[0].shape[0] )

    outlier_indices = np.where( rel_err[0] > error_threshold )

    return indices[0][outlier_indices]

def GetReflectances(ReflectancesFile, KernelsFile,
                    CloudMaskFile,
                    Weigth, InitRow, InitCol,
                    rows, cols, ProcessSnow):

    ReflScaleFactor = 10000.0
    NumberOfBands = 5
    NumberOfParameters = 3

    # BBDR matrix size -> rows x columns x NumberOfBands x 1
    Reflectance = np.zeros((rows, cols, NumberOfBands, 1), np.float32)
    ReflectanceSD = np.zeros((rows, cols, NumberOfBands), np.float32)
    NumberOfSamples = np.zeros((rows, cols), np.float32)

    # Open dataset
    dataset = gdal.Open(ReflectancesFile, GA_ReadOnly)
    band_files = dataset.GetFileList()[1::]

    #SnowMask = dataset.GetRasterBand( ( NumberOfBands * 2 ) + 1 ).ReadAsArray(InitCol,InitRow, cols ,rows)
    # Mask based on ProcessSnow
    #if ProcessSnow == 1:
    #    SnowMask = np.where( SnowMask == 1, 1, 0)
    #else:
    #    SnowMask = np.where( SnowMask == 1, 0, 1)

    # Cloud mask
    d_cloud_mask = gdal.Open(CloudMaskFile)
    CloudMask = d_cloud_mask.ReadAsArray(InitCol,InitRow, cols ,rows).astype(np.bool)

    # Band 5 is 7, no processing bands 5 and 6
    refl_SD_noise_estimates = {1:0.004, 2:0.015, 3:0.003, 4:0.004, 5:0.006 }

    for i in range(NumberOfBands):
        dataset = gdal.Open(band_files[i])
        BandData = dataset.GetRasterBand(1).ReadAsArray(InitCol,InitRow, cols ,rows)
        BandData[BandData == -9999] = 0.0
        #BandData = dataset.GetRasterBand(i+1).ReadAsArray(InitCol,InitRow, cols ,rows)
        #Reflectance[:,:,i,0] = ( BandData / ReflScaleFactor ) * SnowMask * CloudMask
        Reflectance[:,:,i,0] = ( BandData / ReflScaleFactor ) * CloudMask

        #BandData = dataset.GetRasterBand(i+NumberOfBands+1).ReadAsArray(InitCol,InitRow, cols ,rows)
        #ReflectanceSD[:,:,i] = ( BandData / ReflScaleFactor ) * SnowMask * CloudMask
        ReflectanceSD[:,:,i] = refl_SD_noise_estimates[i+1] * CloudMask

    dataset = None

    #-------------------------------------
    # Build C -- coveriance matrix for obs
    #-------------------------------------
    # C in a symetric matrix form of NumberOfBands * NumberOfBands
    C = np.zeros((rows, cols, NumberOfBands, NumberOfBands), np.float32)
    Cinv = np.zeros((rows, cols, NumberOfBands, NumberOfBands), np.float32)

    for i in range(NumberOfBands):
        C[:,:,i,i] = ReflectanceSD[:,:,i] * ReflectanceSD[:,:,i]

    # Create matrices: M, V and E
    # M = Kernels^T C^-1 Kernels
    # V = Kernels^T C^-1 Reflectance
    M = np.zeros((rows, cols, NumberOfBands*NumberOfParameters, NumberOfBands*NumberOfParameters), np.float32)
    V = np.zeros((rows, cols, NumberOfBands*NumberOfParameters), np.float32)

    # Get Kernels
    Kernels = GetKernels(KernelsFile, InitRow, InitCol, rows, cols) 

    for j in range(0,cols):
        for i in range(0,rows):
            if np.where( (Reflectance[i,j,:]>0.0) & (Reflectance[i,j,:]<=1.0) )[0].shape[0] >= 5 and \
               np.where( (ReflectanceSD[i,j,:]>=0.0) & (ReflectanceSD[i,j,:]<=1.0) )[0].shape[0] >= 5 :

                Cinv[i,j,:,:] = np.matrix(C[i,j,:,:]).I
                M[i,j,:,:] = np.matrix(Kernels[i,j,:,:]).T * np.matrix(Cinv[i,j,:,:]) * \
                                       np.matrix(Kernels[i,j,:,:])
                # Multiply only using lead diagonal of Cinv, additionally transpose
                # the result to store V as a 1 x 9 vector
                V[i,j,:] = (np.matrix(Kernels[i,j,:,:]).T * np.diagflat(np.diagonal(Cinv[i,j,:,:])) * Reflectance[i,j,:,:]).T
                NumberOfSamples[i,j] = Weigth

    return ReturnGetReflectances(Reflectance, ReflectanceSD, M, V, NumberOfSamples, Kernels[:,:,0,0:3] )


class ReturnGetReflectances(object):
    def __init__(self, Reflectance, std, M, V, NumberOfSamples, Kernels ):
        self.Reflectance = Reflectance
        self.std = std
        self.M = M
        self.V = V
        self.NumberOfSamples = NumberOfSamples
        self.Kernels = Kernels

def GetKernels(KernelsFile, InitRow, InitCol, rows, cols):

    NumberOfBands = 5
    NumberOfKernels = 3

    Kernels = np.zeros((rows, cols, NumberOfBands, NumberOfKernels * NumberOfBands,), np.float32)
    # Isotropic kernels = 1
    for i in range(NumberOfBands):
         Kernels[:,:,i,i*3] = 1.

    for i in range(1,NumberOfKernels):
        dataset = gdal.Open(KernelsFile, GA_ReadOnly)
        BandData = dataset.GetRasterBand(i).ReadAsArray(InitCol, InitRow, cols ,rows)
        for k in range(NumberOfBands):
            Kernels[:,:,k,(k*3)+i] = BandData

    dataset = BandData = None

    return Kernels


def GetPrior( PriorDataDir, strYear, strDoY, InitRow, InitCol, rows, cols,
              type = 'MCD43', PriorScaleFactor = 10.0 ):

    NumberOfBands = 5
    NumberOfParameters = 3

    # Create matrices
    Prior = np.zeros((rows,cols,NumberOfBands * NumberOfParameters), np.float32)
    PriorVariance = np.zeros((rows,cols,NumberOfBands * NumberOfParameters), np.float32)
    Mask = np.zeros((rows,cols), np.int8)
    NSamples = np.zeros((rows,cols), np.int8)

    C = np.zeros((rows, cols, NumberOfBands * NumberOfParameters, NumberOfBands * NumberOfParameters), np.float32)
    Cinv = np.zeros((rows, cols, NumberOfBands * NumberOfParameters, NumberOfBands * NumberOfParameters), np.float32)
    CinvF = np.zeros((rows, cols, NumberOfBands * NumberOfParameters), np.float32) # Matrix to store C^-1 * Fpr

    if type == 'MCD43':
        if int(strDoY) == 1:
            strDoY = "361"
        else:
            strDoY = '%03d' % ( (int(strDoY) - 8) )

        PriorFile = 'MCD43A1.Prior.%s.img' % ( strDoY )
        PriorFile = os.path.join( PriorDataDir , PriorFile )
        print("Opening prior", PriorFile, "with scaling factor", PriorScaleFactor)
    else:
        # Previous DoY
        Year = int( strYear )
        if int(strDoY) == 1:
            Year = Year - 1
            strYear = str ( Year )
            strDoY = "361"
        else:
            strDoY = '%03d' % ( (int(strDoY) - 8) )

        PriorFile = 'BRDF_Parameters.%s%s.tif' % ( strYear, strDoY )
        PriorFile = os.path.join( PriorDataDir , PriorFile )
        # Check if file exist 
        if glob.glob( PriorFile ) :
            print("Opening prior", PriorFile, "with scaling factor", PriorScaleFactor)
        else:
            return ReturnPriorData( Cinv, CinvF, Prior, PriorVariance, Mask, NSamples )

    dataset = gdal.Open(PriorFile, GA_ReadOnly)
    for i in range(NumberOfBands * NumberOfParameters):
        BandData = dataset.GetRasterBand(i+1).ReadAsArray(InitCol,InitRow, cols ,rows)
        Prior[:,:,i] = BandData

        BandData = dataset.GetRasterBand(i + (NumberOfBands * NumberOfParameters) + 1).ReadAsArray(InitCol,InitRow, cols ,rows)
        # Could be the case that the covariance is 0 or very small but there are samples, make variance = 1
        PriorVariance[:,:,i] = np.where(BandData[:,:] <= 1.0e-8, 1.0, BandData[:,:])
        C[:,:,i,i] = PriorVariance[:,:,i] * PriorScaleFactor
        C[:,:,i,i] = np.where(C[:,:,i,i] > 1.0, 1.0, C[:,:,i,i])

    # Number of Samples 
    nsamples_band = (NumberOfBands*NumberOfParameters*2) + 1
    NSamples[:,:] = dataset.GetRasterBand(nsamples_band).ReadAsArray(InitCol,InitRow, cols ,rows)

    BandData = dataset = None

    # Calculate C inverse
    for j in range(0,cols):
        for i in range(0,rows):
            # Check that al least the isotropic parameters have values
            if np.where( (Prior[i,j,[0,3,6]]>0.0) & (Prior[i,j,[0,3,6]]<=1.0) )[0].shape[0] == 3:

                try:
                    Cinv[i,j,:,:] = np.matrix(C[i,j,:,:]).I
                except np.linalg.LinAlgError:
                    indices = np.where(PriorVariance[i,j,:] == 0.0)[0]
                    PriorVariance[i,j,indices] = 1.0
                    C[i,j,indices,indices] = 1.0
                    Cinv[i,j,:,:] = np.matrix(C[i,j,:,:]).I

                for k in range(NumberOfBands * NumberOfParameters):
                    CinvF[i,j,k] = Cinv[i,j,k,k] * Prior[i,j,k]

                Mask[i,j] = 1

    return ReturnPriorData( Cinv, CinvF, Prior, PriorVariance, Mask, NSamples )

class ReturnPriorData(object):
    def __init__( self, M, V, Parameters, ParametersVariance, Mask, NSamples ):
        self.M = M
        self.V = V
        self.Parameters = Parameters
        self.ParametersVariance = ParametersVariance
        self.Mask = Mask
        self.NSamples = NSamples

def IsLeapYear(Year):
    if Year % 4 == 0:
        if Year % 100 == 0:
            if Year % 400 == 0:
                return True
            else:
                return False
        else:
            return True
    else:
        return False

def GetFileList( strInitDoY, DataDir, Year ):

    # Standard deviation and mean of the normal distribution to create the weighting factors    
    SD = 3.5 #7
    mu = 8.0 #16
    f = NormalDistribution(mu, SD)

    TimeWindow = 16 #32
    FileList = []
    InitDoY = int(strInitDoY)

    # The DoY will be at the center of the 32-day time period
    for Day in range(int(InitDoY-(TimeWindow/2)), \
                     int(InitDoY+(TimeWindow/2))):
        # Some Doy could be lower than 1
        # Assign the right DoY taking into account if Year is a leap-year
        if Day < 1 and IsLeapYear(Year):
            strYear = str(Year - 1)
            strDay = '%03d' % ( 366 + Day )
        elif Day < 1 and not(IsLeapYear(Year)):
            strYear = str(Year - 1)
            strDay = '%03d' % ( 365 + Day )
        # Some DoY could be greater than 365/366
        # Assign the right DoY taking into account if Year is a leap-year
        elif Day == 366 and IsLeapYear(Year):
            strYear = str(Year)
            strDay = '%03d' % ( Day )
        elif Day >= 366 and not(IsLeapYear(Year)):
            strYear = str(Year + 1)
            strDay = '%03d' % ( Day - 365 )
        elif Day >= 367:
            if IsLeapYear(Year):
                strDay = '%03d' % ( Day - 366 )
            else:
                strDay = '%03d' % ( Day - 365 )
            strYear = str(Year + 1)
        else:
            strDay = '%03d' % ( Day )
            strYear = str(Year)

        # Surface reflectance files
        Files = glob.glob( os.path.join( DataDir , strYear,
                           'VRTs' , 'SIAC09GA.%s%s.????.h??v??.vrt' % ( strYear , strDay ) ) )
        for DailyFile in Files:
            FileList.append(DailyFile)

    FileList.sort()

    Year = np.zeros((len(FileList)), np.int16)
    DoY = np.zeros((len(FileList)), np.int16)
    JulianDoY = np.zeros((len(FileList)), np.int16)
    Weigth = np.zeros((len(FileList)), np.float32)

    i = 0
    for j, File in enumerate( FileList ):
        # Get Year and DoY from filename
        YearOfObservation = os.path.basename(File).split('.')[1][0:4]
        DoYOfObservation = os.path.basename(File).split('.')[1][4:7]

        Year[j] = YearOfObservation
        DoY[j] = DoYOfObservation
        JulianDoY[j] = date2num(DT.datetime.strptime(YearOfObservation + DoYOfObservation, "%Y%j"))

        # Observations from the same DoY must have the same weight
        if j == 0:
            DoYWeighting = j
            Weigth[j] = f(DoYWeighting)*((TimeWindow/2)-1)
        else:
            if DoY[j] ==  DoY[j-1]:
                DoYWeighting = i
                Weigth[j] = f(DoYWeighting)*((TimeWindow/2)-1)
            else:
                i += 1
                DoYWeighting = i 
                Weigth[j] = f(DoYWeighting)*((TimeWindow/2)-1)

    # Normalized Weigth
    Weigth = Weigth / Weigth.max()

    return FileList, Year, DoY, JulianDoY, Weigth

def GetDimensions(File):
    dataset = gdal.Open( File, GA_ReadOnly )

    rows, cols = dataset.RasterYSize, dataset.RasterXSize
    dataset = None

    return rows, cols

def Compute_BRDF( LineToProcess, InitCol, NumberOfCols, mask ):
    abs_i = int(LineToProcess)
    abs_j = InitCol

    i = int(LineToProcess - InitRow)

    for j in range( InitCol, NumberOfCols ):
        if mask[abs_i, abs_j] == False:
            abs_j += 1
            continue

        # Outlier detection
        NIR_profile = NIR[i,j,:]
        MinNumberOfObs = 6

        # Perform outlier detection only if there are enough observations
        if NIR_profile.nonzero()[0].shape[0] >= MinNumberOfObs:
            outlier_indices = GetOutlierIndices( NIR_profile, K_profile[i,j] )

            if outlier_indices.shape[0] > 0:
                # Eliminate the observation from the matrices
                M_profile[i,j,:,:, outlier_indices ] = 0.0
                V_profile[i,j,:, outlier_indices]  = 0.0
                tmpNumberOfSamples[i,j, outlier_indices] = 0
                NIR[i,j, outlier_indices ] = 0.0

        NIR_profile = NIR[i,j,:]

        # Assess observations after outlier detection
        valid_obs = AssessObservations( NIR_profile, K_profile[i,j] )
        if valid_obs == False:
            obs_scaler = 0.01
            NumberOfSamples[ abs_i , abs_j ] = 5
        else:
            obs_scaler = 1.0
            NumberOfSamples[ abs_i , abs_j ] = np.sum( tmpNumberOfSamples[i,j,:] )

        # Create accumulators
        M = np.sum( M_profile[i,j,:,:,:], axis = 2 ) * obs_scaler
        V = np.sum( V_profile[i,j,:,:], axis = 1 ) * obs_scaler

        #if (Prior.Mask[i,j] == 1 and NumberOfSamples[ abs_i , abs_j ] >= MinNumberOfObs ):
        if NumberOfSamples[ abs_i , abs_j ] >= MinNumberOfObs:

            # BRDF model inversion without accumulators
            # -----------------------------------------
            # If there is data in the NIR
            indices = np.where( NIR[i,j,:] > 0.0 )
            K = np.ones( [ 3, indices[0].shape[0] ] )
            K[1,:] = K_profile[ i, j, indices, 1 ][0]
            K[2,:] = K_profile[ i, j, indices, 2 ][0]

            K = np.matrix( K )
            #M_ = K * K.T
            #MI = M_.I
            R = NIR[i, j, indices]
            #V_ = K * R.T
            #P_ = ( MI * V_ ).T

            # Get RMSE
            #FWD = P_ * K
            #d = FWD - R
            #e = np.dot( d[0], d[0].T )
            #rmse_ = np.sqrt( e / NumberOfSamples[i,j] )

            #------------------------------------------
            scaler = 1.0
            if PriorPreviousDoY.Mask[i,j] == 1 and \
               PriorPreviousDoY.NSamples[i,j] >= MinNumberOfObs and \
               NumberOfSamples[ abs_i , abs_j ] >= MinNumberOfObs :

                # Enough observations for DoY
                # Enough obs from previous DoY
                #     MELODIES weighted prior from previous DoY
                M_inversion = M + PriorPreviousDoY_Scaled.M[i,j,:,:]
                V_inversion = V + PriorPreviousDoY_Scaled.V[i,j,:]

            elif PriorPreviousDoY.Mask[i,j] == 1 and \
                 PriorPreviousDoY.NSamples[i,j] < MinNumberOfObs and \
                 NumberOfSamples[ abs_i , abs_j ] >= MinNumberOfObs :

                # Enough obs for DoY
                # Not enough obs from previous DoY
                #     MODIS weighted prior for DoY
                M_inversion = M + PriorScaled.M[i,j,:,:]
                V_inversion = V + PriorScaled.V[i,j,:]

            elif PriorPreviousDoY.Mask[i,j] == 1 and \
                 PriorPreviousDoY.NSamples[i,j] >= MinNumberOfObs and \
                 NumberOfSamples[ abs_i , abs_j ] < MinNumberOfObs :

                # Not enough obs for DoY
                # Enough obs from previous DoY
                #     MELODIES non weighted prior from previous DoY
                M_inversion = M + PriorPreviousDoY.M[i,j,:,:]
                V_inversion = V + PriorPreviousDoY.V[i,j,:]

            else:
                # Not enough obs for DoY
                # Not enough obs from previous DoY
                #     MODIS non weighted prior for DoY
                scaler = 1.0 / 20.0
                ##M_inversion = M + Prior.M[i,j,:,:]
                ##V_inversion = V + Prior.V[i,j,:]
                M_inversion = M 
                V_inversion = V

            #(P, rho_residuals, rank, svals) = lstsq(M_inversion, V_inversion)
            P = nnls(M_inversion, V_inversion)
            # TO DO: check under what conditions the BRDF params are negative
            # If negative, use Prior parameters
            #indices = np.where(P < 0.0)
            #P[indices] = Prior.Parameters[i,j,indices]

            # If the prior has some parameters eq 0 is safer to do P for those indices 0
            ##indices = np.where(Prior.Parameters[i,j,:] == 0)
            #P[indices] = 0.0
            ##P[0][indices] = 0.0

            # If the parameters are negative, make them 0.0
            # P[ np.where( P < 0.0 ) ] = 0.0
            #P[0][indices] = 0.0

            #Parameters[ abs_i, abs_j,:,:] = P.reshape(NumberOfBands,NumberOfParameters)
            Parameters[ abs_i, abs_j,:,:] = P[0].reshape( NumberOfBands, NumberOfParameters )

            # Get  NIR RMSE
            #FWD = P[3:6] * K
            FWD = P[0][3:6] * K
            d = FWD - R
            e = np.dot( d[0], d[0].T )
            rmse = np.sqrt( e / NumberOfSamples[ abs_i , abs_j ] )
            GoodnessOfFit[ abs_i, abs_j ] = rmse

            try:
                ParametersVariance[abs_i, abs_j,:,:] = np.diagonal(np.matrix(M_inversion * scaler).I).reshape(NumberOfBands,NumberOfParameters)
            except np.linalg.LinAlgError:
                # Set all off-diagonal elements of M_inversion to 0
                #M_inversion_tmp = np.zeros((M_inversion.shape), np.float32)
                #for k in range(0,9):
                #   M_inversion_tmp[k,k] = M_inversion[k,k]

                #   ParametersVariance[i,j,:,:] = np.diagonal(np.matrix(M_inversion_tmp).I).reshape(NumberOfBands,NumberOfParameters)
                #   M_inversion_tmp = None
                ##Parameters[abs_i,abs_j,:,:] = Prior.Parameters[i,j,:].reshape(NumberOfBands,NumberOfParameters)
                ##ParametersVariance[abs_i,abs_j,:,:] = Prior.ParametersVariance[i,j,:].reshape(NumberOfBands,NumberOfParameters)
                pass

        else:
            pass

            ### If there are no samples at all or there are samples but not prior, use only prior data

            ### Create accumulators
            ##M = np.sum( M_profile[i,j,:,:,:], axis = 2 )
            ##V = np.sum( V_profile[i,j,:,:], axis = 1 )
            ##NumberOfSamples[ abs_i , abs_j ] = np.sum( tmpNumberOfSamples[i,j,:] )

            ##M_inversion = Prior.M[i,j,:,:] + PriorPreviousDoY.M[i,j,:,:] + M
            ##V_inversion = Prior.V[i,j,:] + PriorPreviousDoY.V[i,j,:] + V
            ###(P, rho_residuals, rank, svals) = lstsq(M_inversion, V_inversion)
            ##P = nnls( M_inversion, V_inversion )

            ###Parameters[ abs_i, abs_j,:,:] = P.reshape(NumberOfBands,NumberOfParameters)
            ##Parameters[ abs_i, abs_j, :, :] = P[0].reshape( NumberOfBands, NumberOfParameters )
            ##try:
            ##    scaler = 1.0 / 20.0
            ##    ParametersVariance[abs_i, abs_j,:,:] = np.diagonal(np.matrix(M_inversion * scaler ).I).reshape(NumberOfBands,NumberOfParameters)
            ##except np.linalg.LinAlgError:
            ##    ParametersVariance[abs_i,abs_j,:,:] = Prior.ParametersVariance[i,j,:].reshape(NumberOfBands,NumberOfParameters)

        abs_j += 1


#--------------------------------------------------------------------------------#
from IPython import embed

# Read input 
DoY = int( sys.argv[1] )
strDoY = '%03d' % ( DoY )
strYear = sys.argv[2]

HomeDir = os.path.expanduser('~')

# Get configuration
config = ConfigParser.ConfigParser()
config.read( os.path.join( HomeDir, 'Multiply/src/BRDF/.Multiply_BRDF.cfg' ) )

Project = config.get( 'Project', 'Project' )

cDataDir = config.get( 'Directory', 'DataDir' )
cSrcDir = config.get( 'Directory', 'SrcDir' )
cOutputDir = config.get( 'Directory', 'OutputDir' )
cPriorDataDir = config.get( 'Directory', 'PriorDataDir' )

Sensor = config.get( 'Data', 'Sensor' )
Product = config.get( 'Data', 'Product' )
cTile = config.get( 'Data', 'Tile' )
Format = config.get( 'Data', 'OutputFormat' )

NumberOfBands = config.getint( 'Data', 'NumberOfReflBands' )
NumberOfParameters = config.getint( 'Data', 'NumberOfParameters' )

ProcessSnow = config.getint( 'Processing', 'ProcessSnow' )
NumberOfCores = config.getint( 'Processing', 'NumberOfCores' )

DataDir = os.path.join( HomeDir , Project, cDataDir, Sensor, cTile, Product )

OutputDir = os.path.join( HomeDir , Project, cDataDir, Sensor, cOutputDir, cTile )
PriorDataDir = os.path.join( HomeDir , Project, cDataDir, Sensor, cTile, cPriorDataDir )

print(time.strftime("Processing starting at: %d/%m/%Y %H:%M:%S"))
start = time.time()

# Get list of files to process
FileList, Year, DoY, JulianDoY, Weigth = GetFileList( strDoY, DataDir, int( strYear ) )

# From the first file get dimensions
TotalRows, TotalCols = GetDimensions( FileList[0] )

# From the first file get projection information
tmp_d = gdal.Open( FileList[0] )
gt = tmp_d.GetGeoTransform()
proj = tmp_d.GetProjection()

# Create arrays
Parameters = shm.empty( (TotalRows, TotalCols, 
             NumberOfBands, NumberOfParameters) , np.float32 )
ParametersVariance = shm.empty( (TotalRows, TotalCols, 
                     NumberOfBands, NumberOfParameters), np.float32 )
GoodnessOfFit = shm.empty( ( TotalRows, TotalCols ), np.float32 )
NumberOfSamples = shm.empty( (TotalRows, TotalCols ), np.int16 )

# Depending on the processing system, the composite could be created storing ALL
# datasets in RAM, however for prototyping a tile-based processing will be implemented.
# 100 tiles will be the default setting.
#NumberOfTiles = 12
NumberOfTiles = 80

#iRow = 1
#iCol = 691
#eRow = 2401
#eCol = 891
# Input rows and cols must be in 1-based index
# To process the whole dataset eRow = TotalRows+1, eCol = TotalCols+1
iRow = 1
iCol = 1
eRow = 2401
eCol = 2401

rows = (eRow - iRow)
cols = (eCol - iCol)

for Tile in range(1,NumberOfTiles+1):

    InitRow = (Tile - 1) * (rows / NumberOfTiles) + (iRow - 1)
    EndRow = InitRow + (rows / NumberOfTiles)
    print(InitRow, EndRow)
    InitCol = (iCol - 1)
    print("Processing tile", Tile)

    #Create temporal profile of matrices
    # M = Kernels^T C^-1 Kernels
    # V = Kernels^T C^-1 Reflectance
    # E = Reflectance^T C^-1 Reflectance
    NumberOfFiles = len( FileList )
    NumberOfRows = int(rows / NumberOfTiles)
    NumberOfCols = cols
    M_profile = shm.empty( ( NumberOfRows, NumberOfCols, 
                             NumberOfBands * NumberOfParameters, 
                             NumberOfBands * NumberOfParameters, 
                             NumberOfFiles ), np.float32 )

    V_profile = shm.empty( ( NumberOfRows, NumberOfCols, 
                            NumberOfBands * NumberOfParameters, 
                            NumberOfFiles ), np.float32 )

    K_profile = shm.empty( ( NumberOfRows, NumberOfCols, 
                             NumberOfFiles, 3 ), np.float32 )

    tmpNumberOfSamples = np.zeros( ( NumberOfRows, 
                                     NumberOfCols, 
                                     NumberOfFiles ), np.int16 )

    print((rows / NumberOfTiles))
    NIR = np.zeros((int(rows / NumberOfTiles), cols, NumberOfFiles), np.float32)
    NIR_SD = np.zeros((int(rows / NumberOfTiles), cols, NumberOfFiles), np.float32)

    FileNumber = 0
    for ReflectancesFile in FileList:
        StartReadingTime = time.time()

        print(ReflectancesFile, Weigth[FileNumber])

        KernelsFile = ReflectancesFile[0:len(ReflectancesFile)-3] + "kernels.vrt"
        CloudMaskFile = ReflectancesFile[0:len(ReflectancesFile)-3] + "cloud_mask.vrt"

        Reflectance = GetReflectances( ReflectancesFile,
                          KernelsFile, CloudMaskFile,
                          Weigth[FileNumber], InitRow, InitCol, 
                          int( rows / NumberOfTiles ), cols, ProcessSnow )

        M_profile[:,:,:,:,FileNumber] = Reflectance.M * Weigth[FileNumber]
        V_profile[:,:,:,FileNumber] = Reflectance.V * Weigth[FileNumber]
        K_profile[:,:,FileNumber,:] = Reflectance.Kernels

        tmpNumberOfSamples[:,:,FileNumber] = np.where( Reflectance.NumberOfSamples > 0.0 , 1, 0 )

        #NIR[:,:,FileNumber] = Reflectance.Reflectance[:,:,1,0]
        #NIR_SD[:,:,FileNumber] = Reflectance.std[:,:,1]
        # Blue
        NIR[:,:,FileNumber] = Reflectance.Reflectance[:,:,3,0]
        NIR_SD[:,:,FileNumber] = Reflectance.std[:,:,3]

        Reflectance = None

        FileNumber += 1

        EndReadingTime = time.time()
        TimeElapsed = EndReadingTime - StartReadingTime
        print("Reading data time elapsed = ", (TimeElapsed)/60.0 , "minutes")

    print("Get Land mask...")
    LandMaskFile = "MCD43A2.%s.006.BRDF_Albedo_LandWaterType.tif" % cTile
    LandMaskFile = os.path.join(DataDir, LandMaskFile)
    d_land_mask = gdal.Open(LandMaskFile)
    LandMask = d_land_mask.ReadAsArray()
    LandMask = np.where((LandMask >= 1) & (LandMask <=3), True, False)

    print("Get prior data...")
    #PriorScaled = GetPrior(PriorDataDir, strYear, strDoY, InitRow, InitCol, \
    #    (rows / NumberOfTiles), cols, type = "MCD43", PriorScaleFactor = 20.0)

    # BRDF parameters from previous DoY that would be use as prior
    PriorPreviousDoY_Scaled = GetPrior( OutputDir, strYear, strDoY, \
        InitRow, InitCol, int(rows / NumberOfTiles), cols, type = "MELODIES", PriorScaleFactor = 20.0)

    # Get the prior data to be used when there are less 
    # than 3 samples to perform the BRDF model inversion
    PriorPreviousDoY = GetPrior( OutputDir, strYear, strDoY, \
        InitRow, InitCol, int(rows / NumberOfTiles), cols, type = "MELODIES", PriorScaleFactor = 1.0 )

    # Get the prior to be used only as a gap filler
    #Prior = GetPrior(PriorDataDir, strYear, strDoY, InitRow, InitCol, \
    #    (rows / NumberOfTiles), cols, type = "MCD43", PriorScaleFactor = 1.0)

    # Extract min and max values from NIR_Reflectance
    indices = np.mgrid[0:NIR.shape[0], 0:NIR.shape[1]]
    NIR = np.ma.masked_equal(NIR, 0.0)
    indices_max = NIR.argmax(axis=2)
    indices_min = NIR.argmin(axis=2)

    # Store NIR refl range
    NIR_ReflectanceStats = np.zeros((NumberOfRows, NumberOfCols, 4), np.float32)
    NIR_ReflectanceStats[:,:,0] = NIR[indices[0], indices[1], indices_max]
    NIR_ReflectanceStats[:,:,1] = NIR_SD[indices[0], indices[1], indices_max]
    NIR_ReflectanceStats[:,:,2] = NIR[indices[0], indices[1], indices_min]
    NIR_ReflectanceStats[:,:,3] = NIR_SD[indices[0], indices[1], indices_min]

    print("Performing BRDF model inversion...")

    Processes = []
    NumProcesses = NumberOfCores # Number of cores available to do the processing
    LineToProcess = InitRow

    # Run until all the threads are done, and there is no pixels to process
    while Processes or LineToProcess < EndRow:
        # if we aren't using all the processors AND there are lines left to
        # compute, then spawn another thread
        if ( len( Processes ) < NumProcesses ) and LineToProcess < EndRow:

            p = Process( target = Compute_BRDF,
                args = [ LineToProcess, InitCol, NumberOfCols, LandMask ] )

            p.daemon = True
            p.name = str(LineToProcess)
            p.start()
            Processes.append(p)

            LineToProcess += 1

        # in case that we have the maximum number of threads check
        # if any of them are done.
        else:
            for process in Processes:
                if not process.is_alive():
                    Processes.remove(process)
                    if int(float(process.name)) % 100 == 0:
                        print(process.name, 'processed')

print("Writing results to a file...")
format = "GTiff"
driver = gdal.GetDriverByName(format)
fname = 'BRDF_Parameters.%s%s.tif' % ( strYear, strDoY )
fname = os.path.join( OutputDir, fname )

driver_options = ['COMPRESS=DEFLATE',
                  'BIGTIFF=YES',
                  'PREDICTOR=1',
                  'TILED=YES']

new_dataset = driver.Create( fname, \
  TotalCols, TotalRows, ( NumberOfBands * NumberOfParameters * 2 ) + 2, \
  GDT_Float32, options = driver_options )

k = 1
for i in range(NumberOfBands):
    for j in range(NumberOfParameters):
        new_dataset.GetRasterBand(k).WriteArray(Parameters[:,:,i,j])
        new_dataset.GetRasterBand(k+(NumberOfBands*NumberOfParameters)).WriteArray(ParametersVariance[:,:,i,j])
        k += 1

new_dataset.GetRasterBand((NumberOfBands*NumberOfParameters*2) + 1).WriteArray( NumberOfSamples )
new_dataset.GetRasterBand((NumberOfBands*NumberOfParameters*2) + 2).WriteArray( GoodnessOfFit )

new_dataset.SetProjection( proj )
new_dataset.SetGeoTransform( gt )

new_dataset = None

print(time.strftime("Processing finished at: %d/%m/%Y %H:%M:%S"))
end = time.time()
print("Total time elapsed = ", (end - start)/3600.0, "hours =",  (end - start)/60.0 , "minutes")
