
import os
import glob
import numpy as np
import gdal
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib as mpl

tile = 'h17v05'
Year = 2017
DoY = 160

band = "1"
vmax = 0.025
#band = "2"
#vmax = 0.1
NBand = (5*3) + ((int(band) * 3) - 2)
print(NBand)
strBand = ["", "Red", "NIR"]

# Land Mask
DataDir = "/data/MODIS/%s/SIAC09GA" % tile
fname = "MCD43A2.%s.006.BRDF_Albedo_LandWaterType.tif" % tile
fname = os.path.join(DataDir, fname)
LandMask = gdal.Open(fname).ReadAsArray()
LandMask = np.where((LandMask >= 1) & (LandMask <=3), True, False)

# MOLCI43 using *only* MODIS obs
DataDir = "/data/MODIS/processing/%s" % tile
fname = "BRDF_Parameters.%d%03d.tif" % (Year, DoY)
fname = os.path.join(DataDir, fname)
d = gdal.Open(fname)
uncert_MODIS = d.GetRasterBand(NBand).ReadAsArray()
# Apply Land Mask
uncert_MODIS *= LandMask
uncert_MODIS = np.ma.masked_equal(uncert_MODIS, 0)
# Uncert is variance, get sd
uncert_MODIS = np.sqrt(uncert_MODIS)
print("MODIS:", uncert_MODIS.max())

# MOLCI43 using MODIS and OLCI obs
DataDir = "/data/MOLCI/processing/%s" % tile
fname = "MOLCI43.BRDF_Parameters.%d%03d.tif" % (Year, DoY)
fname = os.path.join(DataDir, fname)
d = gdal.Open(fname)
uncert_MOLCI = d.GetRasterBand(NBand).ReadAsArray()
# Apply Land Mask
uncert_MOLCI *= LandMask
uncert_MOLCI = np.ma.masked_equal(uncert_MOLCI, 0)
# Uncert is variance, get sd
uncert_MOLCI = np.sqrt(uncert_MOLCI)
print("MODIS + OLCI:", uncert_MOLCI.max())

# MOLCI must have more samples, get get max val.
max_val = uncert_MODIS.max()

# Plot NSamples for MOLCI43 using MODIS + OLCI
plot_fname = "Uncert/MOLCI43.MODIS_OLCI.%d%03d.uncert_b%s_f0.png" % (Year, DoY, band)
plt.imsave(plot_fname, uncert_MOLCI, vmin = 0.001, vmax = vmax )

# Plot NSamples for MOLCI43 using only MODIS
plot_fname = "Uncert/MOLCI43.MODIS.%d%03d.uncert_b%s_f0.png" % (Year, DoY, band)
plt.imsave(plot_fname, uncert_MODIS, vmin = 0.001, vmax = vmax)

# Plot colorbar only
# Make a figure and axes with dimensions as desired.
fig = plt.figure(figsize=(8, 1))
ax1 = fig.add_axes([0.05, 0.8, 0.9, 0.15])
# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=0.001, vmax=vmax)
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('Uncertainty %s-f0' % strBand[int(band)])

plt.savefig("Uncert/uncert_%s-f0_colorbar_uncert.png" % strBand[int(band)], dpi=300,
            bbox_inches='tight', transparent=True)

#bins = np.linspace(-10, 10, 100)
#plt.hist(uncert_MODIS, bins, alpha=0.5, label='MODIS')
#plt.hist(y, bins, alpha=0.5, label='MODIS + OLCI')
#plt.legend(loc='upper right')
#plt.savefig("Uncert/hist_uncert_%s-f0_colorbar_uncert.png" % strBand[int(band)], dpi=300,
#            bbox_inches='tight', transparent=True)

