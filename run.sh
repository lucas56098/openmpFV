#!/bin/bash

# hydro variables ------
# file naming
export HYDRO_FOLDER="RT_100" 
export HYDRO_TXT_ADDON="RT"

# mesh options
export HYDRO_N_ROW=50
export HYDRO_CARTESIAN=1    # true or false
export HYDRO_IS_1D=1        # true or false
export HYDRO_STRUCTURE=0    # true or false
export HYDRO_LLOYD_STEPS=0 

# boundary cond
export HYDRO_BOUNDARY_COND=1   # -1 : reflective, 1 : zero gradient
export HYDRO_IS_REPEATING=0     # true or false

# simulation specific
export HYDRO_TOTAL_SIM_TIME=1.5
export HYDRO_CFL=0.1
export HYDRO_TOTAL_SNAPSHOTS=10
export HYDRO_SIM_ORDER=2            # 1st or 2nd
export HYDRO_NUM_THREADS=8         # omp thread number

# initial conditions
# 0 : ShockTube
# 1 : KH 
# 2 : RT 
# 3 : quadshock1 
# 4 : quadshock2 
# 5 : const_flow
export HYDRO_IC_VALUE=0

# ------- Optional: -build to compile hydro before running it -------
BUILD=false
if [ "$1" == "-b" ]; then
    BUILD=true
fi

if $BUILD; then
    # set build folder
    if [ ! -d "build" ]; then
        mkdir build
    fi

    # build the program
    cd build || exit 1
    cmake ../src
    cmake --build . --config Release --target all --
    cd .. || exit 1
    echo "Build completed!"
fi

# ------- Run the program -------
if ! $BUILD; then
    if [ $? -eq 0 ]; then
        echo "Starting program..."
        ./build/hydro
    else
        echo "Build failed!"
        exit 1
    fi
fi
