
import os
import numpy as np
import gdal
import matplotlib.pyplot as plt
from scipy import stats

def hplot(x,y,ax,bar=True,log=True,image=0,
          thresh = 10, xlim=[0,1], ylim=[0,1],
          bins=[128,128]):

    #xyrange = [x.min(),x.max()],[y.min(),y.max()]
    xyrange = xlim,ylim
    min,max = np.min([x.min(),y.min()]),np.max([x.max(),y.max()])

    #xyrange = [0,max],[0,max]
    hh, locx, locy = np.histogram2d( x, y, range=xyrange,
                                    bins=bins)

    hh[hh<thresh] = np.nan
    if log:
      hh = np.log(hh)
    image += np.flipud(hh.T)

    if bar:
      im = ax.imshow(image,cmap='viridis',extent=np.array(xyrange).flatten(),\
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
            print(band, band_number)
            b = g.GetRasterBand(band_number)
            kernels_weights[band, param, :, :] = b.ReadAsArray()
            del(b)

    nsamples_band_number = (NumberOfBands*NumberOfParameters*2) + 1
    mask = g.GetRasterBand(nsamples_band_number).ReadAsArray()
    mask = np.where(mask >= 7, True, False)

    return kernels_weights, mask

def process_mcd43(fname_a1, fname_a2):
    """
    Returns the kernel weights from MCD43, as well as the associated mask
    """

    NumberOfBands = 5
    bands = [1, 2, 3, 4, 7]
    NumberOfParameters = 3

    tmplt_a1 = \
        'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Parameters_Band%d'
    tmplt_a2 = \
        'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Band_Quality_Band%d'
    kernels_weights = np.zeros((NumberOfBands, NumberOfParameters, 2400, 2400))

    mask = np.zeros((NumberOfBands,2400, 2400))
    for i, band in enumerate(bands):
        print(tmplt_a1 % ( fname_a1, band))
        print(tmplt_a2 % ( fname_a2, band))
        g = gdal.Open(tmplt_a1 % ( fname_a1, band))
        tmp_data = g.ReadAsArray()
        tmp_data[tmp_data == 32767] = 0
        kernels_weights[i, :, :, :] = tmp_data * 0.001

        g = gdal.Open(tmplt_a2 % ( fname_a2, band))
        mask[i, :, :] = (g.ReadAsArray() <= 1)
    mask = np.all(mask, axis = 0)

    return kernels_weights, mask

def compare(mcd43, mcd43_mask, molci43, molci43_mask, DoY,
            plot_fname = "MCD43_OLCI43_comparison.png",
            output_dir = "/data/MODIS/processing/h17v05"):
    """
    Plot scatter plots of MCD43 vs MOLCI43
    """
    NumberOfBands = 5
    NumberOfParameters = 3
    #bands = [1, 2, 3, 4, 7]

    band_index = [2, 3, 0, 1, 4]
    band_wv  = [469, 555, 645, 859, 2130]

    fig, axs = plt.subplots(nrows = NumberOfParameters,
                            ncols = NumberOfBands,
                            sharex=True, sharey=True)

    fig.text(0.5, 0.04, 'MCD43 2017 DoY %d' % DoY, ha='center')
    fig.text(0.04, 0.5, 'SYN43 2017 DoY %d' % DoY, va='center', rotation='vertical')

    axs = axs.flatten()

    # Mask
    M = np.where((mcd43_mask == False) | \
                 (molci43_mask == False), False, True)

    k = 0
    for param in range(NumberOfParameters):
        #for band in range(NumberOfBands):
        for i, band in enumerate(band_index):

            mcd43_m = np.where((mcd43[band, param] > 0.0 ) &
                               (mcd43[band, param] < 1.0 ), True, False)

            molci43_m = np.where((molci43[band, param] > 0.0) &
                                 (molci43[band, param] < 1.0), True, False)

            m = M * mcd43_m * molci43_m

            slope, intercept, r_value, p_value, std_err = stats.linregress(
                                                              mcd43[band, param][m],
                                                              molci43[band, param][m])

            hplot(mcd43[band, param][m], molci43[band, param][m], axs[k])
            axs[k].set_title("%dnm $f_%d$ $r^2 = %.02f$" % (band_wv[i], param, r_value ** 2),
                             {'fontsize': 6 })

            plt.yticks([0.0, 0.5, 1.0])
            axs[k].tick_params(labelsize = 7)

            k += 1

    plt.savefig(os.path.join(output_dir, plot_fname),
                dpi=300, bbox_inches="tight")

