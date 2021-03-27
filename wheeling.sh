#!/bin/bash
set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w /io/wheelhouse/
    fi
}

# Downbload boost
wget https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.gz

# Fix boost's detection issue of the python include dir 
ln -fs  /opt/python/cp37-cp37m/include/python3.7m /opt/python/cp37-cp37m/include/python3.7

# Specify python version 

versions=( /opt/python/cp37-cp37m/bin /opt/python/cp38-cp38/bin)
for PYBIN in "${versions[@]}"; do
    export PYBIN
    export BOOST_ROOT="${PYBIN}/../boost/"

    # build boost 
    tar -xzf boost_1_67_0.tar.gz
    cd boost_1_67_0 
    ./bootstrap.sh --with-libraries=python,serialization,iostreams,system,regex --with-python="${PYBIN}/python" --with-python-root="${PYBIN}/../"
    ./b2 install --prefix="${BOOST_ROOT}" -j 20
    cd .. 
    
    # install pip packages
    "${PYBIN}/pip" install cmake-build-extension numpy

    # used in the setup.py
    export PYINC=`echo  ${PYBIN}/../include/python*/ | xargs -n 1 | tail -n 1`
    export LD_LIBRARY_PATH=/rdkit/lib:$BOOST_ROOT/lib

    # copy setup.py for /io
    cp /io/setup.py .

    # build wheel
    "${PYBIN}/python" setup.py bdist_wheel -d /io/wheelhouse_pre/

    # We need to run 'bdist_wheel' a second time to include all libs that were build (I am no sure why)
    "${PYBIN}/python" setup.py bdist_wheel -d /io/wheelhouse_pre/

    # find and include missing libs
    repair_wheel /io/wheelhouse_pre/*.whl

    # clean up
    rm -r /io/wheelhouse_pre/
    rm -r boost_1_67_0
done



