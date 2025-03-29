cmake_minimum_required(VERSION 3.10)
project(hydro_bsc_project VERSION 0.1.0 LANGUAGES C CXX)

# Try to find g++-14 on macOS (Homebrew GCC installation)
find_program(GCC_COMPILER NAMES gcc-14)
find_program(Gpp_COMPILER NAMES g++-14)

# Check if g++-14 is found
if(GCC_COMPILER)
    set(CMAKE_CXX_COMPILER ${Gpp_COMPILER})
    set(CMAKE_C_COMPILER ${GCC_COMPILER})
    message(STATUS "Using GCC compiler (g++-14): ${Gpp_COMPILER}")
else()
    # If g++-14 is not found, use system default and inform user
    message(STATUS "Using default system compiler.")
endif()

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find and enable OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found, linking OpenMP to hydro")
endif()

add_executable(hydro 
    
    # main program
    main.cpp 
    Mesh.tpp
    Solver.tpp
    DG_Solver.tpp

    # different cell types
    cell_types/Point.cpp 
    cell_types/Cell.cpp 
    cell_types/Euler_Cell.cpp
    cell_types/DG_Q_Cell.cpp

    # further utilities
    utilities/Functions.cpp
    utilities/InitialCond.cpp
    
    # imported voronoi_mesh_project
    vmp/Halfplane.cpp 
    vmp/VoronoiCell.cpp 
    vmp/VoronoiMesh.cpp
    )

# Set the name of the compiled program to "hydro"
set_target_properties(hydro PROPERTIES OUTPUT_NAME "hydro")

# Link OpenMP to hydro
if(OpenMP_CXX_FOUND)
    target_link_libraries(hydro PUBLIC OpenMP::OpenMP_CXX)
endif()

# Alternative: Force `-fopenmp` manually (if needed)
target_compile_options(hydro PRIVATE -fopenmp)
target_link_options(hydro PRIVATE -fopenmp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)