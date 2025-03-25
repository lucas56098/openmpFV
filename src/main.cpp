#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <sstream>
#include <sys/resource.h>
#include <omp.h>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "cell_types/DG_Q_Cell.h"
#include "cell_types/Euler_Cell.h"
#include "Mesh.h"
#include "Solver.h"
#include "DG_Solver.h"
#include "vmp/VoronoiMesh.h"
#include "utilities/Functions.h"
#include "utilities/InitialCond.h"
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;


// MAIN :  -------------------------------------------------------------------------------------------------------
int main () {

    // SPECIFICATIONS -------------------------------------
    // files
    string folder_name = getenv("HYDRO_FOLDER") ? getenv("HYDRO_FOLDER") : "KH_lowres";
    string txt_addon = getenv("HYDRO_TXT_ADDON") ? getenv("HYDRO_TXT_ADDON") : "KH";

    // grid
    int N_row = getenv("HYDRO_N_ROW") ? atoi(getenv("HYDRO_N_ROW")) : 10;
    bool cartesian = getenv("HYDRO_CARTESIAN") ? atoi(getenv("HYDRO_CARTESIAN")) : true;
    bool is_1D = getenv("HYDRO_IS_1D") ? atoi(getenv("HYDRO_IS_1D")) : true;
    bool structure = getenv("HYDRO_STRUCTURE") ? atoi(getenv("HYDRO_STRUCTURE")) : false;
    int lloyd_steps = getenv("HYDRO_LLOYD_STEPS") ? atoi(getenv("HYDRO_LLOYD_STEPS")) : 0;

    // boundary
    int boundary_cond = getenv("HYDRO_BOUNDARY_COND") ? atoi(getenv("HYDRO_BOUNDARY_COND")) : -1; // -1 : reflecting, 1 : zero gradient
    bool is_repeating = getenv("HYDRO_IS_REPEATING") ? atoi(getenv("HYDRO_IS_REPEATING")) : false;

    // sim specifics
    double total_sim_time = getenv("HYDRO_TOTAL_SIM_TIME") ? atof(getenv("HYDRO_TOTAL_SIM_TIME")) : 1.0;
    double CFL = getenv("HYDRO_CFL") ? atof(getenv("HYDRO_CFL")) : 0.1;
    int total_snapshots = getenv("HYDRO_TOTAL_SNAPSHOTS") ? atoi(getenv("HYDRO_TOTAL_SNAPSHOTS")) : 20;
    int sim_order = getenv("HYDRO_SIM_ORDER") ? atoi(getenv("HYDRO_SIM_ORDER")) : 2; // 1st or 2nd order
    int num_threads = getenv("HYDRO_NUM_THREADS") ? atoi(getenv("HYDRO_NUM_THREADS")) : 1;
    omp_set_num_threads(num_threads);

    // GRID GENERATION ------------------------------------
    Mesh<Euler_Cell> grid(1);
    
    // INITIAL CONDITIONS ---------------
    // 0 : ShockTube, 1 : KH, 2: RT, 3: quadshock1, 4: quadshock2, 5: const_flow
    int ic_value = getenv("HYDRO_IC_VALUE") ? atoi(getenv("HYDRO_IC_VALUE")) : 0;
    string ic_file_name = getenv("HYDRO_IC_FILE") ? getenv("HYDRO_IC_FILE") : "not-specified";
    
    if (ic_value != 6) {
        grid.generate_grid(cartesian, is_1D, N_row, lloyd_steps, is_repeating, structure);
    } else {
        grid.load_mesh_from_file(ic_file_name, is_repeating);
    }

    Point g_acc;
    if (ic_value == 0) {
        initialize_euler_shock_tube(grid);
    } else if (ic_value == 1) {
        initialize_kelvin_helmholtz(grid);
    } else if (ic_value == 2) {
        g_acc = Point(0, -0.1);
        initialize_rayleigh_taylor(grid, g_acc);
        grid.initialize_boundary_struct(Point(0, 0), 0.25, 1);
        grid.initialize_boundary_struct(Point(0.75, 0), 0.25, 1);
    } else if (ic_value == 3) {
        initialize_quad_shock(grid);
    } else if (ic_value == 4) {
        initialize_quad_shock2(grid);
    } else if (ic_value == 5) {
        initialize_const_flow(grid, Point(1, 0));
    }

    // SIMULATION -----------------------------------------
    Solver<Euler_Cell> solver(&grid);
    
    double t_sim = 0;
    int counter = 0;
    double delT = grid.dt_CFL_euler(CFL);
    int save_iter = (static_cast<int>(total_sim_time/delT)/total_snapshots/10)*10;

    // start timer
    auto start = chrono::high_resolution_clock::now();

    while (t_sim < total_sim_time) {
        if (counter%save_iter == 0 || t_sim + delT >= total_sim_time) {
            
            // save the mesh
            grid.save_mesh(folder_name, cartesian, N_row, sim_order, boundary_cond, is_repeating, total_sim_time, txt_addon, counter, t_sim);

            // timing
            auto now = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = now - start;
            double eta = (elapsed.count() / (t_sim/total_sim_time)) - elapsed.count();
            cout << counter << " : " << delT << " : " << t_sim << ", Time: [" << format_time(elapsed.count()) << "<" << format_time(eta) << "]" << endl;
        }
        
        // calc new timestepping and do euler step
        delT = grid.dt_CFL_euler(CFL);
        solver.euler(delT, boundary_cond, sim_order, g_acc);

        // update time and counter
        t_sim += delT;
        counter += 1;
    }

    // total runtime
    auto final = chrono::high_resolution_clock::now();
    chrono::duration<double> total_time = final - start;
    cout << "Total time: " << total_time.count() << endl;

    // maximum memory
    long long maxrss = get_maxrss_memory();
    

    cout << "done" << endl;
    return 0;    
}


