#!/bin/bash

SRCDIR=$HOME/src/MODIS/BRDF

for DoY in `seq 67 365`
do
    loadavg=`uptime | awk '{print $11}'`
    thisloadavg=`echo $loadavg| awk -F \. '{print $1}'`

    # check load isnt too high  
    if [ "$thisloadavg" -ge "4" ]; then
        sleep 120
    else
        # Check we don't have motre than 8 instances of the prior processing
        cnt=`ps a | grep BRDF_MODIS_Prior.py | grep -v grep | wc -l`

        if [ $cnt -ge 8 ]; then
            wait
            echo "nohup python BRDF_MODIS_Prior.py $DoY 0"
            nohup python $SRCDIR/BRDF_MODIS_Prior.py $DoY 0 &
        else
            # Launch some more processes
            echo "nohup python BRDF_MODIS_Prior.py $DoY 0"
            nohup python $SRCDIR/BRDF_MODIS_Prior.py $DoY 0 &
        fi
    fi
done
