#!/bin/bash

# hydro variables ------
# file naming
export HYDRO_FOLDER="testfolder" 
export HYDRO_TXT_ADDON="test"

# mesh options
export HYDRO_N_ROW=10       # cells in x direction (for voronoi -> N = N_ROW**2, for cartesian N_x = N_row, N_y = N_row * L_y/L_x)
export HYDRO_CARTESIAN=0    # true or false
export HYDRO_IS_1D=0        # true or false
export HYDRO_STRUCTURE=0    # true or false
export HYDRO_LLOYD_STEPS=0 
export HYDRO_BOX_L_X=1
export HYDRO_BOX_L_Y=1

# boundary cond
export HYDRO_BOUNDARY_COND=-1   # -1 : reflective, 1 : zero gradient, 2: wind tunnel
export HYDRO_IS_REPEATING=0     # true or false

# simulation specific
export HYDRO_TOTAL_SIM_TIME=1
export HYDRO_CFL=0.1
export HYDRO_TOTAL_SNAPSHOTS=5
export HYDRO_SIM_ORDER=2            # 1st or 2nd
export HYDRO_NUM_THREADS=4       # omp thread number

# initial conditions
# 0 : ShockTube
# 1 : KH 
# 2 : RT (on a x[0, 0.5] y[0, 1.5] domain)
# 3 : quadshock1 
# 4 : quadshock2 
# 5 : const_flow
# 6 : load from file (also specify directory -> HYDOR_IC_FILE)
#      ! only FV + Voronoi --- overwrites ALL other specified mesh mesh options --- !
export HYDRO_IC_VALUE=0
export HYDRO_IC_FILE="src/files/file_load/v_n30_FV2_BC1_1_5s_test_2_step0.csv"

# ------- Optional: -build to compile hydro before running it -------
BUILD=false
PROFILING=false

if [ "$1" == "-b" ]; then
    BUILD=true
elif [[ "$1" == "-bt" || "$1" == "-tb" ]]; then
    BUILD=true
    PROFILING=true
fi

if $BUILD; then
    # set build folder
    if [ ! -d "build" ]; then
        mkdir build
    else
        rm -fr build
        mkdir build
    fi

    # build the program
    cd build || exit 1

    if $PROFILING; then
        cmake ../src -DENABLE_PROFILING=ON
    else
        cmake ../src
    fi

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
