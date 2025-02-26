#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <omp.h>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "Mesh.h"
#include "Solver.h"

template <typename CellType>
Solver<CellType>::Solver() {}

template <typename CellType>
Solver<CellType>::Solver(Mesh<CellType>* gridin) {
    grid = gridin;
}

// PRIVATE helper functions -----------------------------------------------------------------------
// function to get the normal vector on a edge definde by two points
template <typename CellType>
Point Solver<CellType>::get_normal_vec(Point a, Point b) {

    double Delta_x = b.x - a.x;
    double Delta_y = b.y - a.y;

    double norm = sqrt(Delta_x*Delta_x + Delta_y*Delta_y);

    return Point(-Delta_y/norm, Delta_x/norm);
}

// function to rotate an euler vector where [1] and [2] are the velocities and [0], [3] stay const.
template <typename CellType>
array<double, 4> Solver<CellType>::rotate2Deuler(array<double, 4> vec, double angle) {

    array<double, 4> rot_vec =  {    
                                    vec[0],
                                    cos(angle) * vec[1] - sin(angle) * vec[2],
                                    sin(angle) * vec[1] + cos(angle) * vec[2],
                                    vec[3]
                                };

    return rot_vec;
}


// SOLVERS (update the Mesh by a single step) -----------------------------------------------------
// SOLVER: Euler equations ------------------------------------------------------------------------
// implementation of shallow water equation 2D voronoi using HLL solver
template <typename CellType>
void Solver<CellType>::euler(double dt, int boundary_cond, int sim_order, Point g) {

    // Ensure correct cell type
    if constexpr (is_same_v<CellType, Euler_Cell> == false) {
        std::cerr << "euler step called on Mesh with wrong cell type, you must use Euler_Cells" << std::endl;
        exit(EXIT_FAILURE);
    }

    int num_cells = grid->cells.size();
    std::vector<std::array<double, 4>> Q_new(num_cells);  // Preallocate to avoid race conditions

    std::vector<std::array<std::array<double, 4>, 2>> gradients;
    
    // Compute gradients in parallel if sim_order == 2
    if (sim_order == 2) {
        gradients.resize(num_cells);
        #pragma omp parallel for
        for (int i = 0; i < num_cells; i++) {
            gradients[i] = calc_euler_gradients(i, boundary_cond);
        }
    }

    // Compute new cell values in parallel
    #pragma omp parallel for
    for (int i = 0; i < num_cells; i++) {

        std::array<double, 4> puvE_i_n = {grid->cells[i].rho, grid->cells[i].u, grid->cells[i].v, grid->cells[i].E};
        std::array<double, 4> total_flux = {0, 0, 0, 0};

        // Sum edge flux * edge_length
        for (int j = 0; j < grid->cells[i].edges.size(); j++) {

            // Get value of the neighboring cell
            std::array<double, 4> puvE_j_n = get_puvE_j(puvE_i_n, i, j, boundary_cond);
            std::array<double, 4> puvE_i_ext = puvE_i_n;
            std::array<double, 4> puvE_j_ext = puvE_j_n;
            
            if (sim_order == 2) {
                // Perform linear extrapolation in space and time
                puvE_i_ext = linear_extrapolate_euler(puvE_i_n, gradients, dt, i, j);
                puvE_j_ext = linear_extrapolate_neighbour_euler(i, j, puvE_j_n, puvE_i_ext, gradients, dt, boundary_cond);
            }

            // Compute flux
            std::array<double, 4> F_ij = hll_solver_euler_2D(puvE_i_ext, puvE_j_ext, i, j);

            // Get edge length
            double l_i_j = grid->cells[i].edges[j].length;

            // Accumulate total flux * length
            total_flux[0] += F_ij[0] * l_i_j;
            total_flux[1] += F_ij[1] * l_i_j;
            total_flux[2] += F_ij[2] * l_i_j;
            total_flux[3] += F_ij[3] * l_i_j;
        }

        // Compute new values
        std::array<double, 4> puvE_i_np1;
        double A = grid->cells[i].volume;
        puvE_i_np1[0] = puvE_i_n[0] - (dt / A) * total_flux[0];
        //if (puvE_i_np1[0] < 0.0001) {
        //    puvE_i_np1[0] = 0.0001;
        //}
        puvE_i_np1[1] = (puvE_i_n[0] * puvE_i_n[1] - (dt / A) * total_flux[1] + dt * puvE_i_n[0] * g.x) / puvE_i_np1[0];
        puvE_i_np1[2] = (puvE_i_n[0] * puvE_i_n[2] - (dt / A) * total_flux[2] + dt * puvE_i_n[0] * g.y) / puvE_i_np1[0];
        puvE_i_np1[3] = puvE_i_n[3] - (dt / A) * total_flux[3] + dt * (puvE_i_n[1] * g.x + puvE_i_n[2] * g.y);

        // Store new values (preallocated array prevents race conditions)
        Q_new[i] = puvE_i_np1;
    }

    // Apply updates to grid in parallel
    #pragma omp parallel for
    for (int i = 0; i < num_cells; i++) {
        grid->cells[i].rho = Q_new[i][0];
        grid->cells[i].u = Q_new[i][1];
        grid->cells[i].v = Q_new[i][2];
        grid->cells[i].E = Q_new[i][3];
    }
}



// RIEMANN SOLVERS --------------------------------------------------------------------------------
// HLL Solver for Euler equations 2D unstructured
template <typename CellType>
array<double, 4> Solver<CellType>::hll_solver_euler_2D(array<double, 4> puvE_i_n, array<double, 4> puvE_j_n, int i, int j) {

    // rotate puvE_i, puvE_j into n-frame
    Point n = get_normal_vec(grid->cells[i].edges[j].a, grid->cells[i].edges[j].b);
    double theta = atan2(n.y, n.x);
    array<double, 4> puvE_l = rotate2Deuler(puvE_i_n, -theta);
    array<double, 4> puvE_r = rotate2Deuler(puvE_j_n, -theta);

    // using puvE_l, puvE_r (the rotated ones) now calc F_HLLn using HLL solver
    // calc f_l and f_r
    array<double, 4> f_l = get_flux_f_euler(puvE_l);
    array<double, 4> f_r = get_flux_f_euler(puvE_r);

    // calculate wave speeds
    double gamma = grid->cells[0].get_gamma();
    double SL = min(puvE_l[1] - sqrt((gamma * get_P_ideal_gas(puvE_l))/(puvE_l[0])), puvE_r[1] - sqrt((gamma * get_P_ideal_gas(puvE_r))/(puvE_r[0])));
    double SR = max(puvE_l[1] + sqrt((gamma * get_P_ideal_gas(puvE_l))/(puvE_l[0])), puvE_r[1] + sqrt((gamma * get_P_ideal_gas(puvE_r))/(puvE_r[0])));

    // HLL solver for F
    array<double, 4> F_HLL_n;
    if (SL >= 0) {
        F_HLL_n = f_l;
    } else if (SL < 0 && SR > 0) {
        F_HLL_n[0] = (SR * f_l[0] - SL * f_r[0] + SL * SR * (puvE_r[0] - puvE_l[0])) / (SR - SL);
        F_HLL_n[1] = (SR * f_l[1] - SL * f_r[1] + SL * SR * (puvE_r[0] * puvE_r[1] - puvE_l[0] * puvE_l[1])) / (SR - SL);
        F_HLL_n[2] = (SR * f_l[2] - SL * f_r[2] + SL * SR * (puvE_r[0] * puvE_r[2] - puvE_l[0] * puvE_l[2])) / (SR - SL);
        F_HLL_n[3] = (SR * f_l[3] - SL * f_r[3] + SL * SR * (puvE_r[3] - puvE_l[3])) / (SR - SL);
    } else if (SR <= 0) {
        F_HLL_n = f_r;
    }

    // get F_HLL by rotating F_HLLn back into lab frame
    array<double, 4> F_HLL = rotate2Deuler(F_HLL_n, theta);
    return F_HLL;
}


// GRADIENTS AND EXTRAPOLATION --------------------------------------------------------------------
// calculate gradients for euler cell
template<typename CellType>
array<array<double, 4>, 2> Solver<CellType>::calc_euler_gradients(int i, int boundary_cond) {

    double A = grid->cells[i].volume;
    array<double, 4> U_i = puvE_to_U({grid->cells[i].rho, grid->cells[i].u, grid->cells[i].v, grid->cells[i].E});

    array<array<double, 4>, 2> gradientU;
    gradientU[0] = {0.0, 0.0, 0.0, 0.0};
    gradientU[1] = {0.0, 0.0, 0.0, 0.0};

    // go through all edges
    for (int j = 0; j < grid->cells[i].edges.size(); j++) {

        Point n = get_normal_vec(grid->cells[i].edges[j].a, grid->cells[i].edges[j].b);
        double l_i_j = grid->cells[i].edges[j].length;

        // values of other cell
        array <double, 4> U_j = puvE_to_U(get_puvE_j(U_to_puvE(U_i), i, j, boundary_cond));

        // get c_ij and |r_i_j|
        Point f_mid = get_f_mid(i, j);
        Point r_i = grid->cells[i].seed;
        Point mid = Point(f_mid.x - r_i.x, f_mid.y - r_i.y);

        // distance from seed r_i to neighbour seed r_j
        double r_ij = 2 * (mid.x * n.x + mid.y * n.y);
        // vector from midpoint between r_i and r_j to the midpoint of the edge
        Point c_ij = Point(
                            mid.x - 0.5*r_ij*n.x,
                            mid.y - 0.5*r_ij*n.y
                          );

        // add to sum of gradient
        gradientU[0][0] += (l_i_j/A) * (((U_i[0] + U_j[0])/2.0)*n.x + (U_j[0] - U_i[0]) * (c_ij.x/r_ij));
        gradientU[0][1] += (l_i_j/A) * (((U_i[1] + U_j[1])/2.0)*n.x + (U_j[1] - U_i[1]) * (c_ij.x/r_ij));
        gradientU[0][2] += (l_i_j/A) * (((U_i[2] + U_j[2])/2.0)*n.x + (U_j[2] - U_i[2]) * (c_ij.x/r_ij));
        gradientU[0][3] += (l_i_j/A) * (((U_i[3] + U_j[3])/2.0)*n.x + (U_j[3] - U_i[3]) * (c_ij.x/r_ij));

        gradientU[1][0] += (l_i_j/A) * (((U_i[0] + U_j[0])/2.0)*n.y + (U_j[0] - U_i[0]) * (c_ij.y/r_ij));
        gradientU[1][1] += (l_i_j/A) * (((U_i[1] + U_j[1])/2.0)*n.y + (U_j[1] - U_i[1]) * (c_ij.y/r_ij));
        gradientU[1][2] += (l_i_j/A) * (((U_i[2] + U_j[2])/2.0)*n.y + (U_j[2] - U_i[2]) * (c_ij.y/r_ij));
        gradientU[1][3] += (l_i_j/A) * (((U_i[3] + U_j[3])/2.0)*n.y + (U_j[3] - U_i[3]) * (c_ij.y/r_ij));

    }

    // do slope limiting
    array<array<double, 4>, 2> limited_gradientU;
    //limited_gradientU = slope_limit_tvd(gradientU, i, U_i, boundary_cond, 1.7);
    limited_gradientU = slope_limit_maxmin(gradientU, i, U_i, boundary_cond);

    // return gradient which now can be used to extrapolate the values
    return limited_gradientU;
}


// slope limit euler gradient (like in AREPO -> not TVD)
template <typename CellType>
array<array<double, 4>, 2> Solver<CellType>::slope_limit_maxmin(array<array<double, 4>, 2> gradientU, int i, array<double, 4> U_i, int boundary_cond) {

    // slope limiting
    array<double, 4> a_i = {1.0, 1.0, 1.0, 1.0};

    // go through all edges to calculate correct a_i
    for (int j = 0; j < grid->cells[i].edges.size(); j++) {
        // calculate delta_ij
        array<double, 4> delta_ij;
        Point f_mid = get_f_mid(i, j);
        Point centroid = grid->cells[i].centroid;

        for (int k = 0; k < 4; k++) {
            delta_ij[k] = gradientU[0][k] * (f_mid.x - centroid.x) + gradientU[1][k] * (f_mid.y - centroid.y);
        }

        // calculate psi_ij
        array<double, 4> psi_ij;

        for (int k = 0; k < 4; k++) {

            if (delta_ij[k] > 0) {

                double psi_max = U_i[k];

                //loop through all edges to maximise psi
                for (int l = 0; l < grid->cells[i].edges.size(); l++) {
                    // maximise U_j's
                    array<double, 4> U_j = puvE_to_U(get_puvE_j(U_to_puvE(U_i), i, l, boundary_cond));
                    if (U_j[k] > psi_max) {psi_max = U_j[k];}
                }

                psi_ij[k] = (psi_max - U_i[k])/delta_ij[k];

            } else if (delta_ij[k] < 0) {

                double psi_min = U_i[k];

                // loop through all edges to minimize psi
                for (int l = 0; l < grid->cells[i].edges.size(); l++) {

                    // minimize U_j's
                    array<double, 4> U_j = puvE_to_U(get_puvE_j(U_to_puvE(U_i), i, l, boundary_cond));
                    if (U_j[k] < psi_min) {psi_min = U_j[k];}
                }

                psi_ij[k] = (psi_min - U_i[k])/delta_ij[k];

            } else if (delta_ij[k] == 0) {
                psi_ij[k] = 1;
            }
        }

        // update a_i (minimizing a_i = min_{over j}(1, psi_ij))
        for (int k = 0; k < 4; k++) {
            if (psi_ij[k] < a_i[k]) {a_i[k] = psi_ij[k];}
        }
    }

    // apply slope limiters
    gradientU[0][0] = a_i [0] * gradientU[0][0];
    gradientU[1][0] = a_i [0] * gradientU[1][0];
    gradientU[0][1] = a_i [1] * gradientU[0][1];
    gradientU[1][1] = a_i [1] * gradientU[1][1];
    gradientU[0][2] = a_i [2] * gradientU[0][2];
    gradientU[1][2] = a_i [2] * gradientU[1][2];
    gradientU[0][3] = a_i [3] * gradientU[0][3];
    gradientU[1][3] = a_i [3] * gradientU[1][3];

    // return slope limited gradient
    return gradientU;
}


// slope limit euler gradient (like in TESS -> is TVD)
template <typename CellType>
array<array<double, 4>, 2> Solver<CellType>::slope_limit_tvd(array<array<double, 4>, 2> gradientU, int i, array<double, 4> U_i, int boundary_cond, double theta) {

    // slope limiting
    array<double, 4> a_i = {1.0, 1.0, 1.0, 1.0};

    // go through all edges to calculate correct a_i
    for (int j = 0; j < grid->cells[i].edges.size(); j++) {
        // calculate delta_ij
        array<double, 4> delta_ij;
        Point f_mid = get_f_mid(i, j);
        Point centroid = grid->cells[i].centroid;

        for (int k = 0; k < 4; k++) {
            delta_ij[k] = gradientU[0][k] * (f_mid.x - centroid.x) + gradientU[1][k] * (f_mid.y - centroid.y);
        }

        // calculate psi_ij
        array<double, 4> psi_ij;

        for (int k = 0; k < 4; k++) {

            if (delta_ij[k] > 0) {

                array<double, 4> U_j = puvE_to_U(get_puvE_j(U_to_puvE(U_i), i, j, boundary_cond));
                psi_ij[k] = max(theta*(U_j[k] - U_i[k])/delta_ij[k], 0.0);

            } else if (delta_ij[k] < 0) {

                array<double, 4> U_j = puvE_to_U(get_puvE_j(U_to_puvE(U_i), i, j, boundary_cond));
                psi_ij[k] = max(theta*(U_j[k] - U_i[k])/delta_ij[k], 0.0);

            } else if (delta_ij[k] == 0) {
                psi_ij[k] = 1;
            }
        }

        // update a_i (minimizing a_i = min_{over j}(1, psi_ij))
        for (int k = 0; k < 4; k++) {
            if (psi_ij[k] < a_i[k]) {a_i[k] = psi_ij[k];}
        }
    }

    // apply slope limiters
    gradientU[0][0] = a_i [0] * gradientU[0][0];
    gradientU[1][0] = a_i [0] * gradientU[1][0];
    gradientU[0][1] = a_i [1] * gradientU[0][1];
    gradientU[1][1] = a_i [1] * gradientU[1][1];
    gradientU[0][2] = a_i [2] * gradientU[0][2];
    gradientU[1][2] = a_i [2] * gradientU[1][2];
    gradientU[0][3] = a_i [3] * gradientU[0][3];
    gradientU[1][3] = a_i [3] * gradientU[1][3];

    // return slope limited gradient
    return gradientU;
}


// linearly extrapolate cell i states in space and time towards boundary j
template <typename CellType>
array<double, 4> Solver<CellType>::linear_extrapolate_euler(array<double, 4> puvE_i_n, vector<array<array<double, 4>, 2>> &gradients, double dt, int i, int j) {

    // U_i
    double rho = puvE_i_n[0];
    double u = puvE_i_n[1];
    double v = puvE_i_n[2];
    double E = puvE_i_n[3];
    double P = get_P_ideal_gas(puvE_i_n);
    array<double, 4> U_i = puvE_to_U(puvE_i_n);
    array<array<double, 4>, 2> gradient = gradients[i];

    // d
    Point f_mid = get_f_mid(i, j);
    Point centroid = grid->cells[i].centroid;
    Point d = Point(f_mid.x - centroid.x, f_mid.y - centroid.y);

    // U_i_ext
    array<double, 4> U_i_ext = {
                                    U_i[0] + (d.x * gradient[0][0] - (dt/2.0)*(gradient[0][1]))                                                                 + (d.y * gradient[1][0] - (dt/2.0)*(gradient[1][2])),
                                    U_i[1] + (d.x * gradient[0][1] - (dt/2.0)*((-1*u*u*gradient[0][0]) + (2*u*gradient[0][1])))                                 + (d.y * gradient[1][1] - (dt/2.0)*((-1*u*v*gradient[1][0]) + (v*gradient[1][1]) + (u*gradient[1][2]))),
                                    U_i[2] + (d.x * gradient[0][2] - (dt/2.0)*((-1*u*v*gradient[0][0]) + (v*gradient[0][1]) + (u*gradient[0][2])))              + (d.y * gradient[1][2] - (dt/2.0)*((-1*v*v*gradient[1][0]) + (2*v*gradient[1][2]))),
                                    U_i[3] + (d.x * gradient[0][3] - (dt/2.0)*(((E+P)/(rho))*((-1*u*gradient[0][0]) + gradient[0][1]) + (u*gradient[0][3])))    + (d.y * gradient[1][3] - (dt/2.0)*(((E+P)/(rho))*((-1*v*gradient[1][0]) + gradient[1][2]) + (v*gradient[1][3])))
                               };

    return U_to_puvE(U_i_ext);
}


// linearly extrapolate cell i's neighbours states in space and time towards boundary j
template <typename CellType>
array<double, 4> Solver<CellType>::linear_extrapolate_neighbour_euler(int i, int j, array<double, 4> puvE_j_n, array<double, 4> puvE_i_ext, vector<array<array<double, 4>, 2>> &gradients, double dt, int boundary_cond) {

    if (grid->cells[i].edges[j].is_boundary == false) {

        // index_j is the index of the Cell_j in the cell list
        int index_j = grid->cells[i].edges[j].neighbour->index;
        // index_i is the neighbour index of Cell_i inside of Cell_j
        int index_i;

        // get index_i by looping through all edges of Cell_j
        for (int k = 0; k < grid->cells[index_j].edges.size(); k++) {
            if (grid->cells[index_j].edges[k].is_boundary == false) {
                if (grid->cells[index_j].edges[k].neighbour->index == i) {
                    index_i = k;
                }
            }
        }

        return linear_extrapolate_euler(puvE_j_n, gradients, dt, index_j, index_i);
    } else {
        return get_puvE_j(puvE_i_ext, i, j, boundary_cond);
    }

}



// GET functions and converters -------------------------------------------------------------------
// returns puvE_j vector vor given i, j managing boundary handling
template <typename CellType>
array<double, 4> Solver<CellType>::get_puvE_j(array<double, 4> puvE_i, int i, int j, int boundary_cond) {

    // values of other cell (initalized with sink boundaries for boundary cond = 1)
    array<double, 4> puvE_j = {puvE_i[0], boundary_cond * puvE_i[1], boundary_cond * puvE_i[2], puvE_i[3]};
    // if reflective boundary change velocities accordingly
    if (boundary_cond == -1 || (grid->cells[i].seed.x > 0.1 && grid->cells[i].seed.x < 0.9 && grid->cells[i].seed.y > 0.1 && grid->cells[i].seed.y < 0.9)) {
        Point n = get_normal_vec(grid->cells[i].edges[j].a, grid->cells[i].edges[j].b);

        puvE_j[1] = puvE_i[1] - 2*(puvE_i[1]*n.x + puvE_i[2]*n.y)*n.x;
        puvE_j[2] = puvE_i[2] - 2*(puvE_i[1]*n.x + puvE_i[2]*n.y)*n.y;
    }

    // get actual puvE values if its not a boundary
    if (grid->cells[i].edges[j].is_boundary == false) {
        puvE_j[0] = grid->cells[i].edges[j].neighbour->get_rho();
        puvE_j[1] = grid->cells[i].edges[j].neighbour->get_u();
        puvE_j[2] = grid->cells[i].edges[j].neighbour->get_v();
        puvE_j[3] = grid->cells[i].edges[j].neighbour->get_E();
    }

    return puvE_j;

}


// get Pressure for given euler rho, u, v, E
template <typename CellType>
double Solver<CellType>::get_P_ideal_gas(array<double, 4> puvE) {

    double gamma = grid->cells[0].get_gamma();      // assumes gamma is the same for all cells
    double P = (gamma - 1) * (puvE[3] - (0.5 * puvE[0] * (puvE[1]*puvE[1] + puvE[2]*puvE[2]))) ;
    return P;
}


// get F flux (in x-direction) for given cell state
template <typename CellType>
array<double, 4> Solver<CellType>::get_flux_f_euler(array<double, 4> puvE) {

    double P = get_P_ideal_gas(puvE);

    // calc flux
    array<double, 4> flux = {   puvE[0]*puvE[1],
                                puvE[0]*puvE[1]*puvE[1] + P,
                                puvE[0] * puvE[1] * puvE[2],
                                (puvE[3] + P) * puvE[1]
                            };
    
    return flux;

}


// get midpoint of a edge
template <typename CellType>
Point Solver<CellType>::get_f_mid(int i, int j) {
    return Point((grid->cells[i].edges[j].a.x + grid->cells[i].edges[j].b.x)/2.0, (grid->cells[i].edges[j].a.y + grid->cells[i].edges[j].b.y)/2.0);
}

// convert puveE into U vector
template <typename CellType>
array<double, 4> Solver<CellType>::puvE_to_U(array<double, 4> puvE) {

    array<double, 4> U = {puvE[0], puvE[0]*puvE[1], puvE[0]*puvE[2], puvE[3]};\
    return U;
}
    

// convert U into puveE vector
template <typename CellType>    
array<double, 4> Solver<CellType>::U_to_puvE(array<double, 4> U) {

    array<double, 4> puvE = {U[0], U[1]/U[0], U[2]/U[0], U[3]};
    return puvE;
}    

template <typename CellType>
Solver<CellType>::~Solver() {};