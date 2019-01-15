#!/bin/bash

tile=$1
SRCDIR=$HOME/Multiply/src/BRDF/QL
DATADIR="/data/MODIS/$tile/M?D09GA/2017/Filtered"

for DoY in `seq 160 8 176`
do
    printf -v strDoY "%03d" $DoY

    for file in `ls $DATADIR/*2017$strDoY*h17v05.tif $DATADIR/*predRefl*2017$strDoY*.tif` 
    do
        echo $file
        python $SRCDIR/TrueColor_143_RGB.py $file $tile

        dir=`dirname $file`

        convert Scaled_b1.h??v??.tif Scaled_b4.h??v??.tif Scaled_b3.h??v??.tif \
            -quality 85 -modulate 105,125 -sharpen 4 \
            -font AvantGarde-Book -gravity Southeast \
            -pointsize 40 -fill white -draw 'text 10,18 " '2017$strDoY'' \
            -channel RGB -combine $dir/`basename $file tif`143.RGB.jpg

        convert Scaled_b7.h??v??.tif Scaled_b2.h??v??.tif Scaled_b1.h??v??.tif \
            -quality 85 -modulate 105,125 -sharpen 4 \
            -font AvantGarde-Book -gravity Southeast \
            -pointsize 40 -fill white -draw 'text 10,18 " '2017$strDoY'' \
            -channel RGB -combine $dir/`basename $file tif`721.RGB.jpg

        rm Scaled_b?.h??v??.tif
    done

done

# Create animations
#convert -loop 0 -delay 100 -resize %50 *143.RGB.jpg MCD43A2.Prior.143.RGB_1km.gif
#convert -loop 0 -delay 100 *143.RGB.jpg MCD43A2.Prior.143.RGB_500m.gif

