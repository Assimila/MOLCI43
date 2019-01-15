#!/bin/bash

tile=$1
DATADIR=/data/MODIS
year=2017

DATADIR=$DATADIR/$tile/SIAC09GA/$year
OUTPUTDIR=$DATADIR/VRTs
mkdir -p $OUTPUTDIR

for DoY in `seq 152 184`
do
    strDoY=`printf %03d $DoY`
    date=`date -d "${year}-01-01 +$(( ${DoY} - 1 ))days" +%Y%m%d`
    for image in `ls $DATADIR/MODIS_${tile}_${date}_????_band2_sur.tif`
    do
        image_basename=`basename $image | cut -d_ -f1-4`
        hour_min=`basename $image | cut -d_ -f4-4`
        gdalbuildvrt $OUTPUTDIR/SIAC09GA.$year$strDoY.$hour_min.$tile.vrt \
                     $DATADIR/${image_basename}_band?_sur.tif

        gdal_translate -of VRT $DATADIR/${image_basename}_kernels.tif \
                       $OUTPUTDIR/SIAC09GA.$year$strDoY.$hour_min.$tile.kernels.vrt

        gdal_translate -of VRT $DATADIR/${image_basename}_cloud_mask.tif \
                       $OUTPUTDIR/SIAC09GA.$year$strDoY.$hour_min.$tile.cloud_mask.vrt

    done

done
