#!/bin/bash
echo "Downloading and installing Findent"
git clone git@github.austin.utexas.edu:emj923/Findent-3.1.6.git

cd Findent-3.1.6/findent-3.1.6

sudo make install



pwd=$(pwd)
cd ${pwd%/utilities/code_style/Findent-3.1.6/findent-3.1.6}
echo "Convertering fortran code"

for i in *.f08
do
  echo $i
  findent -i8 <$i> ${i%.*}_1.f08
  mv ${i%.*}_1.f08 $i
done

cd utilities/code_style

echo "Cleaning up"
rm -Rf Findent-3.1.6
echo "FINISHED"
