
import os
import glob
import numpy as np
import gdal
import matplotlib.pyplot as plt
import matplotlib as mpl

tile = 'h17v05'
Year = 2017
DoY = 160

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
NBands, NParams = 5, 3
NSamplesBand = ((NBands * NParams) * 2) + 1
d = gdal.Open(fname)
NSamples_MODIS = d.GetRasterBand(NSamplesBand).ReadAsArray()
# Apply Land Mask
NSamples_MODIS *= LandMask
NSamples_MODIS = np.ma.masked_equal(NSamples_MODIS, 0).astype(np.uint8)

# MOLCI43 using MODIS and OLCI obs
DataDir = "/data/MOLCI/processing/%s" % tile
fname = "MOLCI43.BRDF_Parameters.%d%03d.tif" % (Year, DoY)
fname = os.path.join(DataDir, fname)
NBands, NParams = 5, 3
NSamplesBand = ((NBands * NParams) * 2) + 1
d = gdal.Open(fname)
NSamples_MOLCI = d.GetRasterBand(NSamplesBand).ReadAsArray()
# Apply Land Mask
NSamples_MOLCI *= LandMask
NSamples_MOLCI = np.ma.masked_equal(NSamples_MOLCI, 0).astype(np.uint8)

# MOLCI must have more samples, get get max val.
max_val = NSamples_MOLCI.max()
# Plot NSamples for MOLCI43 using MODIS + OLCI
plot_fname = "MOLCI43.MODIS_OLCI.%d%03d.NSamples.png" % (Year, DoY)
plt.imsave(plot_fname, NSamples_MOLCI, vmin = 7, vmax = max_val)

# Plot NSamples for MOLCI43 using only MODIS
plot_fname = "MOLCI43.MODIS.%d%03d.NSamples.png" % (Year, DoY)
plt.imsave(plot_fname, NSamples_MODIS, vmin = 7, vmax = max_val)

# Plot colorbar only
# Make a figure and axes with dimensions as desired.
fig = plt.figure(figsize=(8, 1))
ax1 = fig.add_axes([0.05, 0.8, 0.9, 0.15])
# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=7, vmax=max_val)
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('Number of samples to perform BRDF model inversion')

plt.savefig("colorbar.png", dpi=300,
            bbox_inches='tight', transparent=True)

