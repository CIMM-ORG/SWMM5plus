#!/bin/bash
pwd=$(pwd)
cd ${pwd%/utilities}

for i in *.f08
do
  findent <$i> ${i%.*}_1.f08
  mv ${i%.*}_1.f08 $i
done
