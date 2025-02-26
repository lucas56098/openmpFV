#include <iostream>
#include <vector>
#include <array>
#include <string>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "Mesh.h"

#ifndef Solver_h
#define Solver_h

template <typename CellType>
class Solver {

public:
    Solver<CellType>();
    Solver<CellType>(Mesh<CellType>* gridin);
    ~Solver<CellType>();

    // equation solvers
    void euler(double dt, int boundary_cond = -1, int sim_order = 2, Point g = Point(0, 0));    

private:

    // grid to solve the equations on
    Mesh<CellType>* grid;

    // helper functions
    Point get_normal_vec(Point a, Point b);
    Point get_f_mid(int i, int j);

        // euler related
        array<double, 4> get_puvE_j(array<double, 4> puvE_i, int i, int j, int boundary_cond);
        array<double, 4> get_flux_f_euler(array<double, 4> puvE);
        array<double, 4> rotate2Deuler(array<double, 4> vec, double angle);
        double get_P_ideal_gas(array<double, 4> puvE);
        array<double, 4> puvE_to_U(array<double, 4> puvE);
        array<double, 4> U_to_puvE(array<double, 4> U);
    
    // riemann solvers
    array<double, 4> hll_solver_euler_2D(array<double, 4> puvE_i_n, array<double, 4> puvE_j_n, int i, int j);

    // second order euler fv
    array<array<double, 4>, 2> calc_euler_gradients(int i, int boundary_cond);
    array<array<double, 4>, 2> slope_limit_maxmin(array<array<double, 4>, 2> gradientU, int i, array<double, 4> U_i, int boundary_cond);
    array<array<double, 4>, 2> slope_limit_tvd(array<array<double, 4>, 2> gradientU, int i, array<double, 4> U_i, int boundary_cond, double theta = 1);
    array<double, 4> linear_extrapolate_euler(array<double, 4> puvE_i_n, vector<array<array<double, 4>, 2>> &gradients, double dt, int i, int j);
    array<double, 4> linear_extrapolate_neighbour_euler(int i, int j, array<double, 4> puvE_j_n, array<double, 4> puvE_i_ext, vector<array<array<double, 4>, 2>> &gradients, double dt, int boundary_cond);
    
};

#include "Solver.tpp"

#endif