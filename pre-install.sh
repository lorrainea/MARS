#! /bin/sh

tar -xvf sdsl-lite.tar.gz
cd sdsl-lite
./install.sh "$(pwd)"/libsdsl
mv libsdsl/ ..
