#!/bin/bash -x

SRCDIR=$HOME/Multiply/src/BRDF/QL
DATADIR=/data/MODIS/processing/mosaics

for DoY in `seq 168 168`
#for DoY in `seq 160 8 176`
do
     printf -v strDoY "%03d" $DoY

    # Extract and scale bands from prior 
    for file in `ls -d $DATADIR/2017$strDoY/b??.tif`
    do
        echo $file
        band=`basename $file | cut -c2-3`
        python $SRCDIR/enhanced_QL.py $file $band
    done

    convert Scaled_b1.tif Scaled_b4.tif Scaled_b3.tif \
            -quality 85 -modulate 105,125 -sharpen 4 \
            -font AvantGarde-Book -gravity Southeast \
            -pointsize 40 -fill white -draw 'text 10,18 " '2017$strDoY'' \
            -channel RGB -combine MOLCI43.2017${strDoY}.143.RGB.jpg

    convert Scaled_b7.tif Scaled_b2.tif Scaled_b1.tif \
            -quality 85 -modulate 105,125 -sharpen 4 \
            -font AvantGarde-Book -gravity Southeast \
            -pointsize 40 -fill white -draw 'text 10,18 " '2017$strDoY'' \
            -channel RGB -combine MOLCI43.2017${strDoY}.721.RGB.jpg

    rm Scaled_b?.tif

done

# Create animations
#convert -loop 0 -delay 100 -resize %50 *143.RGB.jpg MCD43A2.Prior.143.RGB_1km.gif
#convert -loop 0 -delay 100 *143.RGB.jpg MCD43A2.Prior.143.RGB_500m.gif

