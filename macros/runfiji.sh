#!/bin/csh

#cd /groups/apig/tomo/Fiji.app/
#./fiji-linux64 -Xmx24g -Xincgc -XX:MaxPermSize=256m -XX:PermSize=256m -XX:NewRatio=5 #-XX:CMSTriggerRatio=50 -XX:+UseCompressedOops -- -port23

cd /groups/bock/bocklab/people/eric/fiji/
./fiji -Xmx24g -Xincgc -XX:MaxPermSize=256m -XX:PermSize=256m -XX:NewRatio=5 -XX:CMSTriggerRatio=50 -XX:+UseCompressedOops -- -port23

