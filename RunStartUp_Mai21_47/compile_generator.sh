cd HiRE_generator-develop
rm -rf build
mkdir build
cd build
cmake ..
make
cp HiRE_parm ../../
