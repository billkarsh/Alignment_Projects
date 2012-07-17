#!/bin/sh

prfx=/groups/apig/tomo/Eric/eb_slide3_20X/Plate1_0/TIF

sed 's|\(\.*\)|'$prfx'/\1|' < $1 > $2

