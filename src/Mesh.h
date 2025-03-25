#include <vector>
#include <string>
#include "Eigen/Dense"
#include "cell_types/Cell.h"
#include "cell_types/Point.h"
#include "cell_types/Euler_Cell.h"
#include "cell_types/DG_Q_Cell.h"
#include "vmp/VoronoiMesh.h"

#ifndef Mesh_h
#define Mesh_h

template <typename CellType>
class Mesh {

public:
    Mesh<CellType>(int N_bfunc = 1);
    ~Mesh<CellType>();

    // cell list
    vector<CellType> cells;
    bool is_cartesian = false;
    double n_horizontal;
    double n_vertical;
    int N_basisfunc;
    
    // generate mesh
    void generate_grid(bool cartesian, bool is_1D, int N_row, int lloyd_iterations = 0, bool repeating = false, bool structure = false);
    
    // calc cfl
    double dt_CFL_euler(double CFL = 0.4);

    // initial conditions
    void initialize_boundary_struct(Point p0, double l_x, double l_y);

    // DG initial conditions
    // intialize general function
    double step_func2D(Point x, double t, Point p0 = Point(0, 0), Point v = Point(0.5, 0.5), double a = 0.3, double b = 0.3);
    double gaussian2D(Point x, double t, Point p0 = Point(0.5, 0.5), Point v = Point(0.5, 0.5), double A = 1, double sigma = 0.25);
    double legendre_basisfunc2D(Point x, int n);
    Point ksi_to_x(Point ksi, int index);
    Point x_to_ksi(Point x, int index);
    void DG_2D_initialize_step_function(Point p0 = Point(0, 0), Point v = Point(0.5, 0.5), double a = 0.3, double b = 0.3);
    void DG_2D_initialize_gaussian_function(Point p0 = Point(0.5, 0.5), double A = 1, double sigma = 0.25);

    // save mesh
    void save_mesh(string folder_name, bool cartesian, int N_row, int sim_order, int boundary_cond, bool is_repeating, double total_sim_time, string addon, int counter, double t_sim);


private:

    // helper functions
    vector<Point> generate_seed_points(int N, bool fixed_random_seed, double min, int max, int rd_seed, bool sort_pts, int sort_precision, int sort_scheme);
    int get_sort_index(Point pt, int sort_grid_size, int sort_scheme);
    void make_cell_boundary_cell(int i);
    int add_struct(vector<Point>* pts, double dist_a, double safety_dist, string structname);
    void do_lloyd_iterations(vector<Point>* pts, int lloyd_iterations);

    // different grid generations
    void generate_uniform_grid2D(Point start, int n_hor, int n_vert, double distx, double disty, bool repeating = false);
    void generate_vmesh2D(vector<Point> pts, bool repeating = false, bool point_insertion = true);
    void generate_uniform_grid1D(Point start, int n, double dist, bool repeating = false);
    void generate_vmesh1D(vector<Point> pts, bool repeating = false);
};

#include "Mesh.tpp"

#endif