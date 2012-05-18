#!/bin/sh

searchterm="Plate1_0/fixed_images"
replaceterm="Plate1_0/TIF"
inwhich="*.xml"

for file in $(grep -l $searchterm $inwhich)
do
echo $file
mkdir -p _prereplace
cp $file _prereplace
sed -i -e "s|$searchterm|$replaceterm|g" $file
done

