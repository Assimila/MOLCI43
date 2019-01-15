
import os
import glob
import numpy as np
import gdal
import matplotlib.pyplot as plt
from scipy import stats

from compare_MOLCI_MCD import process_mcd43

def hplot(x,y,ax,bar=True,log=True,image=0,
          thresh = 10, xlim=[0,1], ylim=[0,1],
          bins=[128,128]):

    #xyrange = [x.min(),x.max()],[y.min(),y.max()]
    xyrange = xlim,ylim
    min,max = np.min([x.min(),y.min()]),np.max([x.max(),y.max()])

    #xyrange = [0,max],[0,max]
    hh, locx, locy = np.histogram2d(x, y, range = xyrange,
                                    bins=bins)

    hh[hh<thresh] = np.nan
    if log:
      hh = np.log(hh)
    image += np.flipud(hh.T)

    if bar:
      im = ax.imshow(image,cmap='viridis', extent =np.array(xyrange).flatten(),\
                     interpolation='nearest')

      ax.plot(xyrange[0],xyrange[1],'--', color = 'purple', lw=1.5)
      #ax.colorbar(fraction=0.04)
    return image

def process_molci(fname_molci):
    """
    Returns the kernel weights from MOLCI43 and associated mask
    """
    NumberOfBands = 5
    NumberOfParameters = 3

    g = gdal.Open(fname_molci)

    kernels_weights = np.zeros((NumberOfBands, NumberOfParameters, 2400, 2400))

    for band in range(NumberOfBands):
        for param in range(NumberOfParameters):
            band_number = ((band) * NumberOfParameters) + param + 1
            #print(band, band_number)
            b = g.GetRasterBand(band_number)
            kernels_weights[band, param, :, :] = b.ReadAsArray()
            del(b)

    nsamples_band_number = (NumberOfBands*NumberOfParameters*2) + 1
    mask = g.GetRasterBand(nsamples_band_number).ReadAsArray()
    mask = np.where(mask >= 7, True, False)

    return kernels_weights, mask

def get_predicted_refl(kernels_fname, kernel_weights, mcd43_mask = None):
    """
    Computes the predicted reflectance for a specific geometry
    """
    kernels_d = gdal.Open(kernels_fname)
    k_vol = kernels_d.GetRasterBand(1).ReadAsArray()
    k_geo = kernels_d.GetRasterBand(2).ReadAsArray()

    KK = np.moveaxis(np.dstack((np.ones_like(np.array(k_vol)),
                               np.array(k_vol), np.array(k_geo))), -2, 0)

    NumberOfBands = 5

    rows, cols = k_vol.shape
    predicted_ref = np.zeros((NumberOfBands, rows, cols), np.float32)

    for band in range(NumberOfBands):
        predicted_ref[band] = (KK.T * kernel_weights[band]).sum(axis=0)
        if mcd43_mask is not None:
            predicted_ref[band] *= mcd43_mask

    return predicted_ref

def save_file(array, fname, LandMask):
    bands, rows, cols = array.shape
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(fname, rows, cols,
                       bands, gdal.GDT_Float32,
                       options=["TILED=YES", "COMPRESS=DEFLATE"])

    for i in range(bands):
        ds.GetRasterBand(i+1).WriteArray(array[i] * LandMask)

    ds.FlushCache()
    ds = None

def compare(mod09, predicted_refl_molci43,
            predicted_refl_mcd43,
            product, year, LandMask, strDoY,
            output_dir = "/data/MODIS/h17v05/MOD09GA/2017/Filtered"):
    """
    Plot scatter plots of MCD43 vs MOLCI43
    """
    NumberOfBands = 5
    band_index = [2, 3, 0, 1, 4]
    band_wv  = [469, 555, 645, 859, 2130]

    fig, axs = plt.subplots(nrows = 2,
                            ncols = NumberOfBands,
                            sharex=True, sharey=True)

    fig.text(0.5, 0.04, 'MOD09 %d DoY %s' % (year, strDoY), ha='center')
    fig.text(0.04, 0.5, 'Predicted reflectance %d DoY %s' % (year, strDoY), 
             va='center', rotation='vertical')

    axs = axs.flatten()

    k = 0
    for i, band in enumerate(band_index):
        if band == 4:
            mod09_b = mod09[band+2] * 0.0001 * LandMask
        else:
            mod09_b = mod09[band] * 0.0001 * LandMask

        mod09_m = np.where((mod09_b > 0.0) &
                           (mod09_b < 1.0), True, False)

        predicted_refl_m = np.where((predicted_refl_molci43[band] > 0.0) &
                                    (predicted_refl_molci43[band] < 1.0) &
                                    (predicted_refl_mcd43[band] > 0.0) &
                                    (predicted_refl_mcd43[band] < 1.0), True, False)

        m = mod09_m * predicted_refl_m

        hplot(mod09_b[m], predicted_refl_molci43[band][m], axs[k])

        slope, intercept, r_value, p_value, std_err = stats.linregress(mod09_b[m],
                                                          predicted_refl_molci43[band][m])    
        print(slope, intercept, r_value, p_value, std_err)

        plt.text(0.7, 0.95, 'SYN43', horizontalalignment='center',
                 fontsize = 7,
                 verticalalignment='top', transform=axs[k].transAxes)

        axs[k].set_title(r"%dnm $r^2$ = %.02f" % (band_wv[i], r_value ** 2),
                         {'fontsize':7 })

        plt.yticks([0.0, 0.5, 1.0])
        axs[k].tick_params(labelsize = 7)

        k += 1

    for i, band in enumerate(band_index):
        if band == 4:
            mod09_b = mod09[band+2] * 0.0001 * LandMask
        else:
            mod09_b = mod09[band] * 0.0001 * LandMask

        mod09_m = np.where((mod09_b > 0.0) &
                           (mod09_b < 1.0), True, False)

        predicted_refl_m = np.where((predicted_refl_mcd43[band] > 0.0) &
                                    (predicted_refl_mcd43[band] < 1.0) &
                                    (predicted_refl_molci43[band] > 0.0) &
                                    (predicted_refl_molci43[band] < 1.0), True, False)

        m = mod09_m * predicted_refl_m * LandMask

        hplot(mod09_b[m], predicted_refl_mcd43[band][m], axs[k])

        slope, intercept, r_value, p_value, std_err = stats.linregress(mod09_b[m],
                                                          predicted_refl_mcd43[band][m])
        print(slope, intercept, r_value, p_value, std_err)

        plt.text(0.7, 0.95, 'MCD43', horizontalalignment='center',
                 fontsize = 7,
                 verticalalignment='top', transform=axs[k].transAxes)

        axs[k].set_title(r'%dnm $r^2$ = %.02f' % (band_wv[i], r_value ** 2),
                         {'fontsize': 7 })

        plt.yticks([0.0, 0.5, 1.0])
        axs[k].tick_params(labelsize = 7)

        k += 1

    plt.savefig(os.path.join(output_dir,
                             "%s_vs_SIAC09.%d%s.png" % (product, year, strDoY)),
                dpi=300, bbox_inches="tight")

tile = 'h17v05'
product = 'MYD09GA'
year = 2017
mDataDir = '/data/MODIS'

tile_mask_fname = "MCD43A2.%s.006.BRDF_Albedo_LandWaterType.tif" % tile
tile_mask_fname = os.path.join(mDataDir, tile, 'SIAC09GA', tile_mask_fname)
d_land_mask = gdal.Open(tile_mask_fname)
LandMask = d_land_mask.ReadAsArray()
LandMask = np.where((LandMask >= 1) & (LandMask <=3), True, False)

for DoY in [160]:
    strDoY = "%03d" % int(DoY)

    # Get MCD43
    print("Getting MCD43...")
    fname_a1 = 'MCD43A1.A%d%s.%s.006.*.hdf' % (year, strDoY, tile)
    fname_a1 = os.path.join(mDataDir, tile, 'MCD43A1', fname_a1)
    fname_a1 = glob.glob(fname_a1)[0]
    print(fname_a1)

    fname_a2 = 'MCD43A2.A%d%s.%s.006.*.hdf' % (year, strDoY, tile)
    fname_a2 = os.path.join(mDataDir, tile, 'MCD43A2', fname_a2)
    fname_a2 = glob.glob(fname_a2)[0]
    print(fname_a2)

    mcd43, mcd43_mask = process_mcd43(fname_a1, fname_a2)


    # Get MOD09
    print("Getting MOD09...")
    mod09_fname = "%s.%d%s.%s.tif" % (product, year, strDoY, tile)
    DataDir = os.path.join(mDataDir, tile, product, str(year), 'Filtered')
    mod09_fname = os.path.join(DataDir, mod09_fname)
    d_mod09 = gdal.Open(mod09_fname)
    mod09 = d_mod09.ReadAsArray()

    kernels_fname = "%s.%d%s.%s.kernels.tif" % (product, year, strDoY, tile)
    kernels_fname = os.path.join(DataDir, kernels_fname)

    # Get MOLCI43
    print("Getting MOLCI43...")
    #molci_fname = 'BRDF_Parameters.%d%s.tif' % (year, strDoY)
    mDataDir = '/data/MOLCI'
    molci_fname = 'MOLCI43.BRDF_Parameters.%d%s.tif' % (year, strDoY)
    molci_fname = os.path.join(mDataDir, 'processing', tile, molci_fname)
    molci43, molci43_mask = process_molci(molci_fname)

    # Get predicted reflectances
    predicted_refl_molci43 = get_predicted_refl(kernels_fname, molci43)
    predicted_refl_mcd43 = get_predicted_refl(kernels_fname, mcd43, mcd43_mask)

    # Save predicted refl
    mDataDir = '/data/MODIS'
    DataDir = os.path.join(mDataDir, tile, product, str(year), 'Filtered/MOLCI43_MODIS_and_OLCI')
    output_fname = os.path.join(DataDir,
                       'MOLCI43.%s.predRefl.%d%s.tif' % (product, year, strDoY))
    save_file(predicted_refl_molci43, output_fname, LandMask)

    output_fname = os.path.join(DataDir,
                       'MCD43.%s.predRefl.%d%s.tif' % (product, year, strDoY))
    save_file(predicted_refl_mcd43, output_fname, LandMask)

    # Plot
    compare(mod09, predicted_refl_molci43,
            predicted_refl_mcd43,
            product, year, LandMask,
            strDoY, output_dir = DataDir)

