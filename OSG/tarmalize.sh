#!/bin/sh
WD=`pwd`
echo "Compiling urlencode. . . "
$CXX urlencode.cpp -o urlencode
cd ..
echo "Compacting input data. . . "
./IC79_diffuse --writeCompact --exitAfterLoading
if [ "$?" -ne 0 ]; then
	cd $WD
	echo "Failed to compact data"
	exit 1
fi
cd $WD
mkdir -p osg_diffuse
cd osg_diffuse
echo "Gathering files. . . "
cp -p ../../IC79_diffuse ./
cp -p ../../compact_data.dat ./
cp -pr ../../dom_eff_fits ./
cp -pr ../../../../NewNuFlux/resources/data ./flux_paramterizations
cp -pr ../urlencode ./
mkdir -p lib
for lib in `ldd IC79_diffuse | sed -n 's|.*=> \(.*\) (.*)|\1|p'`; do
	cp -p $lib ./lib/
done
cd ..
echo "Compressing tarball. . . "
tar cf osg_diffuse.tar osg_diffuse
if [ -f osg_diffuse.tar.xz ]; then
	rm osg_diffuse.tar.xz
fi
xz -9 -M 4096MiB -T 1 -e osg_diffuse.tar
echo "Cleaning up. . . "
rm -rf osg_diffuse
echo "Done"
