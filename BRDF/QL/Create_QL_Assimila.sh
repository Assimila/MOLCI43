#!/bin/bash -x

tile=$1
SRCDIR=$HOME/Multiply/src/BRDF/QL
DATADIR=/data/MOLCI/processing/$tile
#DATADIR=/data/MODIS/processing/$tile

mask=/data/MODIS/$tile/SIAC09GA/MCD43A2.$tile.006.BRDF_Albedo_LandWaterType.tif

#for DoY in `seq 1 8 365`
for DoY in `seq 160 8 176`
#for DoY in `seq 137 8 297`
do
     printf -v strDoY "%03d" $DoY

    # Extract and scale bands from prior 
    #for file in `ls -d $DATADIR/BRDF_Parameters.2017$strDoY.tif`
    for file in `ls -d $DATADIR/MOLCI43.BRDF_Parameters.2017$strDoY.tif`
    do
        echo $file
        tile=`echo $file | cut -d/ -f5`
        python $SRCDIR/Prior_enhanced_QL.py $file $tile $tile
    done

    convert Scaled_b1.h??v??.tif Scaled_b4.h??v??.tif Scaled_b3.h??v??.tif \
            -quality 85 -modulate 105,125 -sharpen 4 \
            -font AvantGarde-Book -gravity Southeast \
            -pointsize 40 -fill white -draw 'text 10,18 " '2017$strDoY'' \
            -channel RGB -combine MOLCI43.MODIS.2017${strDoY}.143.RGB.jpg

    convert Scaled_b7.h??v??.tif Scaled_b2.h??v??.tif Scaled_b1.h??v??.tif \
            -quality 85 -modulate 105,125 -sharpen 4 \
            -font AvantGarde-Book -gravity Southeast \
            -pointsize 40 -fill white -draw 'text 10,18 " '2017$strDoY'' \
            -channel RGB -combine MOLCI43.MODIS.2017${strDoY}.721.RGB.jpg

    rm Scaled_b?.h??v??.tif

done

# Create animations
convert -loop 0 -delay 100 -resize %50 *143.RGB.jpg MOLCI43.MODIS.2017.143.RGB_1km.gif
convert -loop 0 -delay 100 -resize %50 *721.RGB.jpg MOLCI43.MODIS.2017.721.RGB_1km.gif

#convert -loop 0 -delay 100 -resize %50 *143.RGB.jpg MCD43A2.Prior.143.RGB_1km.gif
#convert -loop 0 -delay 100 *143.RGB.jpg MCD43A2.Prior.143.RGB_500m.gif

