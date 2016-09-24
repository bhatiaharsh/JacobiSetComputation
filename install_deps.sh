#! /bin/bash

IN_DIR=$PWD

# --------------------------------------
# install trimesh2

if [ ! -d "trimesh2" ]; then
wget http://gfx.cs.princeton.edu/proj/trimesh2/src/trimesh2-2.12.tar.gz
tar -xzf trimesh2-2.12.tar.gz
fi

cd trimesh2
rm -r lib.*
rm -r bin.*
make -j8

cd $IN_DIR
rm trimesh2-2.12.tar.gz

# --------------------------------------
# install and patch SoS

if [ ! -d "Detri_2.6.a" ]; then
wget http://www.geom.uiuc.edu/software/cglist/GeomDir/Detri_2.6.a.tar.gz
tar -xzf Detri_2.6.a.tar.gz
fi

patch -s -p0 < patch_SOS.txt

cd Detri_2.6.a
mkdir build
cd build

cmake ../
make -j8

cd $IN_DIR
rm Detri_2.6.a.tar.gz