#!/bin/sh

export MRC_TRIM=12

cd lou_example
scr green_raw.xml -dtemp -p_N_

cd temp

cp /groups/apig/tomo/lou_stuff/thmparams.txt .
cp /groups/apig/tomo/lou_stuff/dmeshparams.txt .
cp /groups/apig/tomo/lou_stuff/subfm.sh .
cp /groups/apig/tomo/lou_stuff/sub8.sh .
cp /groups/apig/tomo/lou_stuff/thmdown.sh .
cp /groups/apig/tomo/lou_stuff/thmsame.sh .
cp /groups/apig/tomo/lou_stuff/combine.sh .
cp /groups/apig/tomo/lou_stuff/finish.sh .
cp /groups/apig/tomo/lou_stuff/report.sh .

cp /groups/apig/tomo/lou_stuff/runlsq.sh stack/

./subfm.sh 0 2
