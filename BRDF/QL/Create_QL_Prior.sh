#!/bin/bash -x

SRCDIR=$HOME/src/MODIS/BRDF/QL
#DATADIR=$HOME/MELODIES/data/MODIS
DATADIR=/data/MODIS

tile=$1

for DoY in `seq 1 365`
do
     printf -v strDoY "%03d" $DoY

    # Extract and scale bands from prior 
    for file in `ls -d $DATADIR/$tile/Prior/MCD43A1.Prior.$strDoY.img`
    do
        echo $file
        #tile=`echo $file | cut -d/ -f4`
        python $SRCDIR/Prior_enhanced_QL.py $file $tile
    done

    convert Scaled_b1.h??v??.tif Scaled_b4.h??v??.tif Scaled_b3.h??v??.tif \
            -quality 85 -modulate 105,125 -sharpen 4 \
            -font AvantGarde-Book -gravity Southeast \
            -pointsize 40 -fill white -draw 'text 10,18 " '$strDoY'' \
            -channel RGB -combine MCD43A2.Prior.${strDoY}.143.RGB.jpg

    rm Scaled_b?.h??v??.tif

done

# Create animations
convert -loop 0 -delay 100 -resize %50 *143.RGB.jpg MCD43A2.Prior.143.RGB_1km.gif
convert -loop 0 -delay 100 *143.RGB.jpg MCD43A2.Prior.143.RGB_500m.gif

