#!/bin/sh
mkdir DEScan2
cp NAMESPACE DEScan2/
cp DESCRIPTION DEScan2/
cp NEWS DEScan2/
cp -rv R DEScan2/
cp -rv inst DEScan2/
cp -rv src DEScan2/
#cp -rv testData DEScan2/
cp -rv tests DEScan2/
#cp -rv vignettes DEScan2/
cp -rv man DEScan2/

cd DESCan2
find ./ -name ".DS_Store" -depth -exec rm {} \;
cd testData
echo "removing ..."
rm -rv new_files
cd ../src
rm -v *.*o
rm -v *.rds
echo "...done!"
cd ../..
#R CMD check --no-vignettes DEScan2
R CMD BiocCheck --no-check-vignettes DEScan2
