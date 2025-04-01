#!/bin/bash

# load eigen lib
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar -xzvf eigen-3.4.0.tar.gz

# copy into src
cp -r ./eigen-3.4.0/Eigen ./src/

# cleanup
rm eigen-3.4.0.tar.gz
rm -fr eigen-3.4.0