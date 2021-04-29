yum install -y wget freetype-devel zlib-devel libpng12-devel pixman-devel

# compile cairo (RDKit needs older version than in centOS7 repo)
wget https://www.cairographics.org/releases/cairo-1.10.0.tar.gz --no-check-certificate
tar xvf cairo-*
cd cairo-*
./configure
make -j 20
make install
cd ..
