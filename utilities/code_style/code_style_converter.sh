#!/bin/bash
echo "Downloading and installing Findent"
git clone git@github.austin.utexas.edu:emj923/Findent-3.1.6.git


cd Findent-3.1.6/findent-3.1.6
pwd=$(pwd)

./configure --prefix=${pwd%/utilities/code_style/Findent-3.1.6/findent-3.1.6}
sudo make install


pwd=$(pwd)
cd ${pwd%/utilities/code_style/Findent-3.1.6/findent-3.1.6}


echo "Convertering fortran code"
for i in *.f08
do
  echo $i
  ./bin/findent -i4 <$i> ${i%.*}_1.f08
  mv ${i%.*}_1.f08 $i
done


echo "Cleaning Up"


sudo rm -Rf bin
sudo rm -Rf share


cd utilities/code_style


rm -Rf Findent-3.1.6
echo "FINISHED"
