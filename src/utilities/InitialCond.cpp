#include "InitialCond.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <sys/resource.h>

// function to initalize sods shock tube on the mesh
void initialize_euler_shock_tube(Mesh<Euler_Cell> &grid) {

    vector<Euler_Cell>& cells = grid.cells;

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        if (cells[i].seed.x < 0.5) { // + cells[i].seed.y < 1
            cells[i].rho += 1-1;
            cells[i].u = 0;
            cells[i].v = 0;
            cells[i].E += (1/(cells[i].gamma - 1)) * 1 - 1;
        } else {
            cells[i].rho += 0.125-1;
            cells[i].u = 0;
            cells[i].v = 0;
            cells[i].E += (1/(cells[i].gamma - 1)) * 0.1 - 1;
        }
    }
}

// function to initalize a kelvin helmholtz instability on the mesh
void initialize_kelvin_helmholtz(Mesh<Euler_Cell> &grid) {

    vector<Euler_Cell>& cells = grid.cells;

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        double pi = 3.14159265358979323846;
        double sigma = 0.05/sqrt(2);

        if (cells[i].seed.y > 0.25 && cells[i].seed.y < 0.75) {
            cells[i].rho = 2; //0.5
            cells[i].u = 0.5; //0.3
            cells[i].v = 0.1*sin(4*pi*cells[i].seed.x)*(exp(-((cells[i].seed.y -0.25)*(cells[i].seed.y -0.25))/(2 * sigma * sigma)) + exp(-((cells[i].seed.y -0.75)*(cells[i].seed.y -0.75))/(2 * sigma * sigma)));
            cells[i].E = (2.5/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u);
        } else {
            cells[i].rho = 1; //0.2
            cells[i].u = -0.5; //-0.3
            cells[i].v = 0.1*sin(4*pi*cells[i].seed.x)*(exp(-((cells[i].seed.y -0.25)*(cells[i].seed.y -0.25))/(2 * sigma * sigma)) + exp(-((cells[i].seed.y -0.75)*(cells[i].seed.y -0.75))/(2 * sigma * sigma)));
            cells[i].E = (2.5/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u);
        }
    }
  
}


// function to initalize a rayleigh taylor instability on the mesh
void initialize_rayleigh_taylor(Mesh<Euler_Cell> &grid, Point g) {

    vector<Euler_Cell>& cells = grid.cells;

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        double pi = 3.14159265358979323846;

        cells[i].gamma = 1.4;

        if (cells[i].seed.y > 0.5){// + 0.03*cos(pi*2*((cells[i].seed.x + 0.25)*2))) {
            cells[i].rho = 2;
            cells[i].u = 0;
            cells[i].v = 0.0025*(1-cos(4*pi * (cells[i].seed.x - 0.5)))*(1 - cos(2*pi*cells[i].seed.y));
            double P = 2.5 + cells[i].rho * g.y * (cells[i].seed.y - 0.5);
            cells[i].E = (P/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u);
        } else {
            cells[i].rho = 1;
            cells[i].u = 0;
            cells[i].v = 0.0025*(1-cos(4*pi * (cells[i].seed.x - 0.5)))*(1 - cos(2*pi*cells[i].seed.y));
            double P = 2.5 + cells[i].rho * g.y * (cells[i].seed.y - 0.5);
            cells[i].E = (P/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u);
        }
    }
  
}


// function to initialize a quad shock
void initialize_quad_shock(Mesh<Euler_Cell> &grid) {

    vector<Euler_Cell>& cells = grid.cells;

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        if (cells[i].seed.x >= 0.5 && cells[i].seed.y >= 0.5) {
            cells[i].rho = 1.5;
            cells[i].u = 0;
            cells[i].v = 0;
            cells[i].E = (1.5/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        } else if (cells[i].seed.x < 0.5 && cells[i].seed.y > 0.5) {
            cells[i].rho = 0.5323;
            cells[i].u = 1.206;
            cells[i].v = 0;
            cells[i].E = (0.3/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        } else if (cells[i].seed.x < 0.5 && cells[i].seed.y < 0.5) {
            cells[i].rho = 0.138;
            cells[i].u = 1.206;
            cells[i].v = 1.206;
            cells[i].E = (0.029/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        } else if (cells[i].seed.x > 0.5 && cells[i].seed.y < 0.5) {
            cells[i].rho = 0.5323;
            cells[i].u = 0;
            cells[i].v = 1.206;
            cells[i].E = (0.3/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        }
    }

}


// function to initialize another quad shock
void initialize_quad_shock2(Mesh<Euler_Cell> &grid) {

    vector<Euler_Cell>& cells = grid.cells;

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        cells[i].gamma = 7.0/5.0;

        if (cells[i].seed.x >= 0.5 && cells[i].seed.y >= 0.5) {
            cells[i].rho = 1;
            cells[i].u = 0.75;
            cells[i].v = -0.5;
            cells[i].E = (1/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        } else if (cells[i].seed.x < 0.5 && cells[i].seed.y > 0.5) {
            cells[i].rho = 2;
            cells[i].u = 0.75;
            cells[i].v = 0.5;
            cells[i].E = (1/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        } else if (cells[i].seed.x < 0.5 && cells[i].seed.y < 0.5) {
            cells[i].rho = 1;
            cells[i].u = -0.75;
            cells[i].v = 0.5;
            cells[i].E = (1/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        } else if (cells[i].seed.x > 0.5 && cells[i].seed.y < 0.5) {
            cells[i].rho = 3;
            cells[i].u = -0.75;
            cells[i].v = -0.5;
            cells[i].E = (1/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        }
    }

}


// function to initialize a constant flow
void initialize_const_flow(Mesh<Euler_Cell> &grid, Point v) {
    
    vector<Euler_Cell>& cells = grid.cells;
    
    for (int i = 0; i<cells.size(); i++) {
        cells[i].u += v.x;
        cells[i].v += v.y;
    }
}