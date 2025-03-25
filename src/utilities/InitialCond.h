#include "../cell_types/Point.h"
#include "../Mesh.h"
#include <vector>
#include <string>
using namespace std;

#ifndef InitialCond_h
#define InitialCond_h

void initialize_euler_shock_tube(Mesh<Euler_Cell> &grid);
void initialize_kelvin_helmholtz(Mesh<Euler_Cell>&grid);
void initialize_rayleigh_taylor(Mesh<Euler_Cell> &grid, Point g);
void initialize_quad_shock(Mesh<Euler_Cell> &grid);
void initialize_quad_shock2(Mesh<Euler_Cell> &grid);
void initialize_const_flow(Mesh<Euler_Cell> &grid, Point v);

#endif