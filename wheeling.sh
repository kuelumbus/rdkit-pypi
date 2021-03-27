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

# Download boost 
wget https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.gz

# Fix: Boost has an issue detecting the include directory of python3.7m (because of the "m")
ln -fs  /opt/python/cp37-cp37m/include/python3.7m /opt/python/cp37-cp37m/include/python3.7

# Python versions
versions=(/opt/python/cp37-cp37m/bin /opt/python/cp38-cp38/bin /opt/python/cp39-cp39/bin)
for PYBIN in "${versions[@]}"; do
    export PYBIN
    export BOOST_ROOT="${PYBIN}/../boost/"

    # Build boost
    tar -xzf boost_1_67_0.tar.gz
    cd boost_1_67_0 
    ./bootstrap.sh --with-libraries=python,serialization,iostreams,system,regex --with-python="${PYBIN}/python" --with-python-root="${PYBIN}/../"
    ./b2 install --prefix="${BOOST_ROOT}" -j 20
    cd .. 
    
    # Install packages necessary for the build
    "${PYBIN}/pip" install cmake-build-extension numpy

    # Used later in setup.py
    export PYINC=`echo  ${PYBIN}/../include/python*/ | xargs -n 1 | tail -n 1`
    export LD_LIBRARY_PATH=/rdkit/lib:$BOOST_ROOT/lib

    # Copy setup.py from /io
    cp /io/setup.py .

    # Build wheel
    "${PYBIN}/python" setup.py bdist_wheel -d /io/wheelhouse_pre/

    # We need to run 'bdist_wheel' a second time to include all libs that were built (I am not sure why)
    "${PYBIN}/python" setup.py bdist_wheel -d /io/wheelhouse_pre/

    # Find and include missing libs
    repair_wheel /io/wheelhouse_pre/*.whl

    # Clean up
    rm -r /io/wheelhouse_pre/
    rm -r boost_1_67_0
done



