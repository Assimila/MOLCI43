
import os
import glob
from compare_MOLCI_MCD import process_mcd43, process_molci, compare

year = 2017
tile = 'h17v05'

mDataDir = '/data/MODIS'

#for DoY in range(160, 176+1, 8):
for DoY in range(160, 161):

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

    # Get MOLCI43
    mDataDir = '/data/MOLCI'
    DataDir = os.path.join(mDataDir, 'processing', tile)
    #fname_molci = 'BRDF_Parameters.%d%s.tif' % (year, strDoY)
    fname_molci = 'MOLCI43.BRDF_Parameters.%d%s.tif' % (year, strDoY)
    fname_molci = os.path.join(DataDir, fname_molci)
    print(fname_molci)
    
    molci43, molci43_mask = process_molci(fname_molci)

    # Plot
    compare(mcd43, mcd43_mask, molci43, molci43_mask, DoY = DoY,
            plot_fname = 'MCD43_SYN43_MODIS_and_OLCI_comparison.%d%s.%s.png' % (year, strDoY, tile),
            output_dir = '/data/MOLCI/processing/h17v05')
