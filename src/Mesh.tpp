#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <string>
#include "Mesh.h"
#include "Eigen/Dense"
#include "cell_types/Cell.h"
#include "cell_types/Point.h"
#include "cell_types/Euler_Cell.h"
#include "vmp/VoronoiMesh.h"
#include "utilities/Functions.h"
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <type_traits>
#include <random>
#include <sstream>

template <typename CellType>
Mesh<CellType>::Mesh() {
    N_basisfunc = 0;
}

template <typename CellType>
Mesh<CellType>::Mesh(int N_bfunc) {
    N_basisfunc = N_bfunc;
}


// PRIVATE Helper functions for Point Generation --------------------------------------------------
// Function out of vmp main.cpp
// RANDOM POINTS: function to get a sort index
template <typename CellType>
int Mesh<CellType>::get_sort_index(Point pt, int sort_grid_size, int sort_scheme) {

    double nr = static_cast<double>(sort_grid_size);

    int index;

    // sort by x-y modulo grid
    if (sort_scheme == 1) {

        if (static_cast<int>(nr * nr * pt.y)%2==0) {

            index = (sort_grid_size - static_cast<int>(pt.x * nr)) + static_cast<int>(nr * nr * pt.y);

        } else {

            index = static_cast<int>(pt.x * nr) + static_cast<int>(nr * nr * pt.y);

        }

    // sort radially outward
    } else if (sort_scheme == 2) {

        index = static_cast<int>((nr*nr)*sqrt((pt.x-0.5)*(pt.x-0.5) + (pt.y-0.5)*(pt.y-0.5)));

    // sort radially inward
    } else if (sort_scheme == 3) {

        index = static_cast<int>((nr*nr)*(sqrt(0.5) - sqrt((pt.x-0.5)*(pt.x-0.5) + (pt.y-0.5)*(pt.y-0.5))));

    // all other numbers -> do not sort
    } else {

        index = 1;

    }

    return index;
}


// Function out of vmp main.cpp
// RANDOM POINTS: generates seed points to use for mesh generation
/*
template <typename CellType>
vector<Point> Mesh<CellType>::generate_seed_points(int N, bool fixed_random_seed, double min, int max, int rd_seed, bool sort_pts, int sort_precision, int sort_scheme) {
    vector<Point> points;

    unsigned int random_seed;
    default_random_engine eng;

    // set either fixed or changing random seed
    if (fixed_random_seed) {
        random_seed = rd_seed;
    } else {
        random_device rd;
        random_seed = rd();

    }

    // define uniform random distribution
    eng = default_random_engine(random_seed);
    uniform_real_distribution<double> distr(min, max);
    
    /*
    // optional point density for KH-instability
    uniform_real_distribution<int> distra(min, max);
    normal_distribution<double> distr1(0.7, 0.07);
    normal_distribution<double> distr2(0.3, 0.07);


    // generate random coordinates for Points
    for (int i = 0; i < N; ++i) {
        double x = distr(eng);

        double a;
        if (distra(eng) > 0) {
            a = distr1(eng);
        } else {
            a = distr2(eng);
        }

        double y = std::max(std::min(a, 0.9999), 0.0001) + 0.0009999999 * distr(eng);

        points.push_back(Point(x, y));
    }
    */

    /*
    // optional point density for quad shock
    uniform_real_distribution<int> distra(1, 2);
    for (int i = 0; i < N; i++) {
        double x;
        double y;
        double scalea = 0.2*distr(eng) + 0.6;
        double scaleb = 0.2*distr(eng) + 0.6;
        if (distra(eng) % 10 != 0) {
            x = scalea*distr(eng) + 0;
            y = scaleb*distr(eng) + 0;
        } else {
            x = distr(eng);
            y = distr(eng);
        }

        points.push_back(Point(x,y));
    }*/

    /*
    // optional point density for circle
    normal_distribution<double> distr2(0.5, 0.08);

    // generate random coordinates for Points
    for (int i = 0; i < N; ++i) {
        double x = distr(eng);

        double a = distr2(eng);


        double y = std::max(std::min(a, 0.9999), 0.0001) + 0.0009999999 * distr(eng);

        points.push_back(Point(x, y));
    }*/
    
    /*

   
    // generate random coordinates for Points
    for (int i = 0; i < N; ++i) {
        double x = distr(eng);
        double y = distr(eng);

        points.push_back(Point(x, y));
    }
    

    // if this is true the points will be sorted
    if (sort_pts) {
        vector<int> indices;
        vector<int> sort_indices;

        // get sort indices
        for (int i = 0; i < points.size(); i++) {
            indices.push_back(get_sort_index(points[i], sort_precision, sort_scheme));
            sort_indices.push_back(i);
        }

        // combine data into pairs
        vector<pair<int, int> > combined;
        for (int i = 0; i < indices.size(); ++i) {
            combined.push_back(make_pair(indices[i], sort_indices[i]));
        }    

        // sort combined data by sort indices
        sort(combined.begin(), combined.end());

        // get sorted_pts
        vector<Point> sorted_pts;
        for (int i = 0; i < combined.size(); i++) {
            sorted_pts.push_back(points[combined[i].second]);
        }

        return sorted_pts;
    }

    return points;
}
*/


// RANDOM POINTS: generates seed points to use for mesh generation
template <typename CellType>
vector<Point> Mesh<CellType>::generate_seed_points(int N, bool fixed_random_seed, Point min, Point max, int rd_seed, bool sort_pts, int sort_precision, int sort_scheme) {
    vector<Point> points;

    unsigned int random_seed;
    default_random_engine eng;

    // set either fixed or changing random seed
    if (fixed_random_seed) {
        //cout << "specify random_seed: ";
        //cin >> random_seed;
        random_seed = rd_seed;
    } else {
        random_device rd;
        random_seed = rd();

    }

    // define uniform random distribution
    eng = default_random_engine(random_seed);
    uniform_real_distribution<double> distrx(min.x, max.x);
    uniform_real_distribution<double> distry(min.y, max.y);

    // generate random coordinates for Points
    for (int i = 0; i < N; ++i) {
        double x = distrx(eng);
        double y = distry(eng);
        points.push_back(Point(x, y));
    }

    // if this is true the points will be sorted
    if (sort_pts) {
        vector<int> indices;
        vector<int> sort_indices;

        // get sort indices
        for (int i = 0; i < points.size(); i++) {
            indices.push_back(get_sort_index(points[i], sort_precision, sort_scheme));
            sort_indices.push_back(i);
        }

        // combine data into pairs
        vector<pair<int, int> > combined;
        for (int i = 0; i < indices.size(); ++i) {
            combined.push_back(make_pair(indices[i], sort_indices[i]));
        }    

        // sort combined data by sort indices
        sort(combined.begin(), combined.end());

        // get sorted_pts
        vector<Point> sorted_pts;
        for (int i = 0; i < combined.size(); i++) {
            sorted_pts.push_back(points[combined[i].second]);
        }

        return sorted_pts;
    }

    return points;
}

// GRID GENERATION: -------------------------------------------------------------------------------
// calls the generate Mesh functions depending on specified options (cartesian, 1D/2D, N_row, optional lloyd preprocessing, repeating boundary conditions)
template <typename CellType>
void Mesh<CellType>::generate_grid(bool cartesian, bool is_1D, int N_row, int lloyd_iterations, bool repeating, bool structure, double L_x, double L_y) {

    if (cartesian) {
        // generate cartesian mesh
        if (is_1D) {
            // do it in 1D
            this->generate_uniform_grid1D(Point(0, 0), N_row, L_x/static_cast<double>(N_row), repeating);
            is_cartesian = true;
        } else {
            // do it in 2D
            this->generate_uniform_grid2D(Point(0, 0), N_row, static_cast<int>(static_cast<double>(N_row)*(L_y/L_x)), L_x/static_cast<double>(N_row), L_x/static_cast<double>(N_row), repeating);
            is_cartesian = true;
        }
    } else {
        // generate voronoi mesh
        if (is_1D) {
            // do it in 1D
            vector<Point> pts = generate_seed_points(N_row, true, Point(0, 0), Point(L_x, L_y), 42, true, 100, 1);
            this->generate_vmesh1D(pts, repeating); // rescaling does not work here (l_x and l_y = 1 by default here)
            is_cartesian = false;
        } else {
            // do it in 2D
            vector<Point> pts = generate_seed_points(N_row * N_row, true, Point(0, 0), Point(L_x, L_y), 42, true, 100, 1);
            if (lloyd_iterations != 0) {do_lloyd_iterations(&pts, lloyd_iterations, L_x, L_y);};
            int nr;
            //int nr2;
            if (structure) {nr = add_struct(&pts, 0.02, 0.05, "struct");}
            //if (structure) {nr2 = add_struct(&pts, 0.000001, 0.006, "struct_europe");}

            this->generate_vmesh2D(pts, repeating, !structure, L_x, L_y);
            is_cartesian = false;

            // counting here is for each structure: 0-nr: innner points, nr-2nr: outer points, 2nr - ... olt pts
            if (structure) {
                for (int i = 0; i < nr; i++) { // if activating second structure here nr needs to be changed to nr2
                    make_cell_boundary_cell(i);
                }
                //for (int i = 2*nr2; i < 2*nr2 + nr; i++) {
                //    make_cell_boundary_cell(i);
                //}
            }

        }

    }
    for (int i = 0; i<this->cells.size(); i++) {
        this->cells[i].index = i;
    }
    cout << "grid generated" << endl;
    


}


// generates a uniform grid with all the neighbour relations and so on
template <typename CellType>
void Mesh<CellType>::generate_uniform_grid2D(Point start, int n_hor, int n_vert, double distx, double disty, bool repeating) {

    n_horizontal = n_hor;
    n_vertical = n_vert;

    for (int b = 0; b<n_vert; b++) {
        for (int a = 0; a<n_hor; a++) {

            // set the correct seed
            Point seedin(start.x + 0.5*distx + a*distx, start.y + 0.5*disty + b*disty);
            
            // start to define the faces
            vector<edge> edgesin;
            edge f0; edge f1; edge f2; edge f3;

            // set a and b for the faces
            f0.a = Point(start.x + a*distx, start.y + b*disty);
            f0.b = Point(start.x + a*distx, start.y + b*disty + disty);
            f1.a = Point(start.x + a*distx, start.y + b*disty + disty);
            f1.b = Point(start.x + a*distx + distx, start.y + b*disty + disty);
            f2.a = Point(start.x + a*distx + distx, start.y + b*disty + disty);
            f2.b = Point(start.x + a*distx + distx, start.y + b*disty);
            f3.a = Point(start.x + a*distx + distx, start.y + b*disty);
            f3.b = Point(start.x + a*distx, start.y + b*disty);

            // set correct boundary flags
            f0.is_boundary = false; f1.is_boundary = false; f2.is_boundary = false; f3.is_boundary = false;
            if (a == 0 && repeating == false) {f0.is_boundary = true;}
            if (a == n_hor -1 && repeating == false) {f2.is_boundary = true;}
            if (b == 0 && repeating == false) {f3.is_boundary = true;}
            if (b == n_vert -1 && repeating == false) {f1.is_boundary = true;}

            if (n_vert == 1) {
                f1.is_boundary = true;
                f3.is_boundary = true;
            }

            // push faces in edges vector
            edgesin.push_back(f0);
            edgesin.push_back(f1);
            edgesin.push_back(f2);
            edgesin.push_back(f3);

            // set edge length to dist
            edgesin[0].length = disty;
            edgesin[1].length = distx;
            edgesin[2].length = disty;
            edgesin[3].length = distx;

            // push back new cell in cells vector
            if constexpr (is_same_v<CellType, DG_Q_Cell> == true) {
                cells.emplace_back(seedin, edgesin, N_basisfunc);
            } else {
                cells.emplace_back(seedin, edgesin);
            }

            // set centroid
            cells[cells.size()-1].centroid = seedin;

        }
    }

    // now that all cells exist we define the neighbour relations
    for (int i = 0; i<cells.size(); i++) {

        cells[i].volume = distx * disty;

        for (int j = 0; j<cells[i].edges.size(); j++) {
            if (cells[i].edges[j].is_boundary == false) {
                if (j == 0) {
                    cells[i].edges[j].neighbour = &cells[i - (i%n_hor) + ((i+n_hor - 1)%n_hor)];
                } else if (j == 1) {
                    cells[i].edges[j].neighbour = &cells[((((i - (i%n_hor))/n_hor)+1)%n_vert)*n_hor + (i%n_hor)];
                } else if (j == 2) {
                    cells[i].edges[j].neighbour = &cells[i - (i%n_hor) + ((i + 1)%n_hor)];
                } else if (j == 3) {
                    cells[i].edges[j].neighbour = &cells[((((i - (i%n_hor))/n_hor)+n_vert - 1)%n_vert)*n_hor + (i%n_hor)];
                }
                
            }
        }
    }
}


// generates vmesh using vmp and converts it into data usable for this mesh type
template <typename CellType>
void Mesh<CellType>::generate_vmesh2D(vector<Point> pts, bool repeating, bool point_insertion, double L_x, double L_y) {

    // preprocessing for repeating boundary conditions
    vector<Point> points_plus_ghost;
    int initial_pts_size = pts.size();
    if (repeating) {
        points_plus_ghost.reserve(pts.size() * 9);

        // put points into (middle/middle) block by shrinking them by a factor of 3
        for (int i = 0; i<pts.size(); i++) {
            points_plus_ghost.emplace_back((pts[i].x/3.0) + 1.0/3.0, (pts[i].y/3.0) + 1.0/3.0);
        }

        // add the same shrinked points again but shifted in all other 8 third blocks (up/middle/down, left/middle/right)
        vector<double> pos_X = {0., 1., 2., 0., 2., 0., 1., 2.};
        vector<double> pos_Y = {0., 0., 0., 1., 1., 2., 2., 2.};
        for (int i = 0; i<8; i++) {
            for (int j = 0; j<pts.size(); j++) {
                points_plus_ghost.emplace_back((pts[j].x/3.0) + pos_X[i] * 1.0/3.0, (pts[j].y/3.0) + pos_Y[i] * 1.0/3.0);
            }
        }
        
        // replace pts with pts + additional ghost cells (eg 8 times the pts all around)
        pts = points_plus_ghost;
    }

    // generate vmesh
    VoronoiMesh vmesh(pts, L_x, L_y);

    if (point_insertion) {
        vmesh.do_point_insertion();
    } else {
        vmesh.construct_mesh();
    }


    // loop through all cells (of initial pts vector) to set everything but neighbour relations
    for (int i = 0; i<initial_pts_size; i++) {

        // set seed and define edge vector
        Point seedin(vmesh.vcells[i].seed.x, vmesh.vcells[i].seed.y);
        vector<edge> edgesin;

        for (int j = 0; j<vmesh.vcells[i].edges.size(); j++) {

            edge f;

            // set edge quantities
            f.a = vmesh.vcells[i].verticies[((vmesh.vcells[i].edges.size()-1) + j)%vmesh.vcells[i].edges.size()];
            f.b = vmesh.vcells[i].verticies[j];
            f.length = sqrt((f.a.x - f.b.x)*(f.a.x - f.b.x) + (f.a.y - f.b.y)*(f.a.y - f.b.y));
            f.is_boundary = (vmesh.vcells[i].edges[j].index1 < 0);

            edgesin.push_back(f);

        }

        // push back new cell in cells vector
        if constexpr (is_same_v<CellType, DG_Q_Cell> == true) {
            cells.emplace_back(seedin, edgesin, N_basisfunc);
        } else {
            cells.emplace_back(seedin, edgesin);
        }

        // set centroid
        cells[i].centroid = vmesh.vcells[i].get_centroid();

    }

    // loop through all cells (of initial pts size) and set neighbour relations
    for (int i = 0; i<initial_pts_size; i++) {

        cells[i].volume = vmesh.vcells[i].get_area();

        for (int j = 0; j<vmesh.vcells[i].edges.size(); j++) {

            // neighbour index as used in vmesh, using modulo here to get neighbour relations
            // correct even for repeating boundaries since index will be multiple of original index
            int neigbour_index;
            neigbour_index = (vmesh.vcells[i].edges[j].index2)%initial_pts_size;

            // exclude boundaries
            if (neigbour_index >= 0) {

                // neighbour as defined before over index now with adress
                cells[i].edges[j].neighbour = &cells[neigbour_index];

            }

        }

    }

    // if repeating boundary conditions rescale the mesh back to normal
    if (repeating) {

        // loop through all cells
        for (int i = 0; i<cells.size(); i++) {

            // redo scaling for seed and volume
            cells[i].seed.x = (cells[i].seed.x - 1.0/3.0) * 3.0;
            cells[i].seed.y = (cells[i].seed.y - 1.0/3.0) * 3.0;
            cells[i].centroid.x = (cells[i].centroid.x - 1.0/3.0) * 3.0;
            cells[i].centroid.y = (cells[i].centroid.y - 1.0/3.0) * 3.0;
            cells[i].volume = cells[i].volume * 3.0 * 3.0;

            // loop through all faces
            for (int j = 0; j<cells[i].edges.size(); j++) {

                edge f = cells[i].edges[j];

                // redo scaling for edge positions
                f.a.x = (f.a.x - 1.0/3.0) * 3.0;
                f.a.y = (f.a.y - 1.0/3.0) * 3.0;
                f.b.x = (f.b.x - 1.0/3.0) * 3.0;
                f.b.y = (f.b.y - 1.0/3.0) * 3.0;
                cells[i].edges[j].a.x = f.a.x;
                cells[i].edges[j].a.y = f.a.y;
                cells[i].edges[j].b.x = f.b.x;
                cells[i].edges[j].b.y = f.b.y;

                // recalculate length with correct scaling
                cells[i].edges[j].length = sqrt((f.a.x - f.b.x)*(f.a.x - f.b.x) + (f.a.y - f.b.y)*(f.a.y - f.b.y));

            }

        }

    }


}


// does the lloyd_iterations to some points
template <typename CellType>
void Mesh<CellType>::do_lloyd_iterations(vector<Point>* pts, int lloyd_iterations, double L_x, double L_y) {
    // preprocessing step to change pts for mesh into pts for approx centroidal vmesh
    if (lloyd_iterations != 0) {
        cout << "start lloyd_iterations" << endl;

        // calculate original mesh
        VoronoiMesh initial_vmesh(*pts, L_x, L_y);
        initial_vmesh.do_point_insertion();
        
        // do multiple iterations of lloyds algorithm
        for (int i = 0; i<lloyd_iterations; i++) {

            // calculate centroids
            vector<Point> centroids;
            centroids.reserve(initial_vmesh.vcells.size());
            for (int i = 0; i<initial_vmesh.vcells.size(); i++) {
                centroids.push_back(initial_vmesh.vcells[i].get_centroid());
            }

            // replace original mesh seeds with centroids and calculate mesh again
            initial_vmesh = VoronoiMesh(centroids, L_x, L_y);
            initial_vmesh.do_point_insertion();

        }

        // after iterations store final calculated seeds in pts
        vector<Point> centroidal_seeds;
        centroidal_seeds.reserve(initial_vmesh.vcells.size());
        for (int i = 0; i<initial_vmesh.vcells.size(); i++) {
            centroidal_seeds.push_back(initial_vmesh.vcells[i].seed);
        }
        *pts = centroidal_seeds;
        cout << "finished lloyd_iterations" << endl;
    }
}


// generates a 1D uniform grid with all the neighbour relations and so on
template <typename CellType>
void Mesh<CellType>::generate_uniform_grid1D(Point start, int n, double dist, bool repeating) {
    generate_uniform_grid2D(start, n, 1, dist, 1, repeating);
}


// generates a 1D voronoi mesh, only works if points are between 0 and 1
template <typename CellType>
void Mesh<CellType>::generate_vmesh1D(vector<Point> pts, bool repeating) {

    for (int i = 0; i < pts.size(); i++) {
        pts[i].y = 0.5;
    }

    vector<Point> sorted_pts = pts;
    sort(sorted_pts.begin(), sorted_pts.end(), [](const Point& a, const Point& b) {
        return a.x < b.x;
    });

    for (int i = 0; i < sorted_pts.size(); i++) {

        // set the correct seed
        Point seedin = sorted_pts[i];

        // start to define the faces
        vector<edge> edgesin;
        edge f0; edge f1; edge f2; edge f3;

        // set a and b for the faces and boundary flags
        f0.is_boundary = false;
        f1.is_boundary = true;
        f2.is_boundary = false;
        f3.is_boundary = true;

        double distl;
        double distr;

        // manually set left part of faces
        if (i == 0) {
            // we are in the leftmost cell
            f0.is_boundary = true;
            f0.a = Point(0, 0);
            f0.b = Point(0, 1);
            f1.a = Point(0, 1);
            f3.b = Point(0,0);
            distl = seedin.x;
        } else {
            // just a random cell
            f0.a = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 0);
            f0.b = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 1);
            f1.a = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 1);
            f3.b = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 0);
            distl = sorted_pts[i].x - sorted_pts[i-1].x;
        }
        
        // manually set rigth part of faces
        if (i == sorted_pts.size()-1) {
            // we are in the rightmost cell
            f2.is_boundary = true;
            f1.b = Point(1, 1);
            f2.a = Point(1, 1);
            f2.b = Point(1,0);
            f3.a = Point(1,0);
            distr = 1 - seedin.x;
        } else {
            // just a random cell
            f1.b = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 1);
            f2.a = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 1);
            f2.b = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 0);
            f3.a = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 0);
            distr = sorted_pts[i+1].x - sorted_pts[i].x;
        }

        if (repeating == true) {
            f0.is_boundary = false;
            f2.is_boundary = false;
        }

        // push back faces
        edgesin.push_back(f0);
        edgesin.push_back(f1);
        edgesin.push_back(f2);
        edgesin.push_back(f3);

        // set correct lengths
        edgesin[0].length = 1;
        edgesin[2].length = 1;
        edgesin[1].length = distl + distr;
        edgesin[3].length = distl + distr;

        // push back new cell in cells vector
        if constexpr (is_same_v<CellType, DG_Q_Cell> == true) {
            cells.emplace_back(seedin, edgesin, N_basisfunc);
        } else {
            cells.emplace_back(seedin, edgesin);
        }

        cells[cells.size()-1].centroid = seedin;

    }

    // now that all cells exist we define the neighbour relations
    for (int i = 0; i < cells.size(); i++) {

        cells[i].volume = cells[i].edges[0].length * cells[i].edges[1].length;

        if (cells[i].edges[0].is_boundary == false) {
            cells[i].edges[0].neighbour = &cells[i - (i%cells.size()) + ((i+cells.size() - 1)%cells.size())];
        }
        if (cells[i].edges[2].is_boundary == false) {
            cells[i].edges[2].neighbour = &cells[i - (i%cells.size()) + ((i + 1)%cells.size())];
        }

    }

}


// Boundary Structures -------------------------------------------------------------------------
// Function to create an internal boundary structure (just a square at the moment)
template <typename CellType>
void Mesh<CellType>::initialize_boundary_struct(Point p0, double l_x, double l_y) {

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        // if cell[i] is inside square make it boundary cell (in principle any other if condition could be built here)
        if (cells[i].seed.x < p0.x + l_x && cells[i].seed.x > p0.x && cells[i].seed.y < p0.y + l_y && cells[i].seed.y > p0.y) { 
            make_cell_boundary_cell(i);
        }
    }
}


// makes cell[i] boundary cell
template <typename CellType>
void Mesh<CellType>::make_cell_boundary_cell(int i) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Euler_Cell> == true) {
        // boundary cells have h = -INFINITY such that they will not be plotted in visualization
        cells[i].rho = -INFINITY;
        cells[i].u = 0;
        cells[i].v = 0;
        cells[i].E = -INFINITY;
    }

    // loop through edges of that boundary cell
    for (int j = 0; j < cells[i].edges.size(); j++) {

        // if already boundary there is nothing to change
        if (cells[i].edges[j].is_boundary == false) {

            // set edge from internal side to boundary
            cells[i].edges[j].is_boundary = true;

            // go through all faces of the edges[j].neighbour were looking at
            for (int k = 0; k < cells[i].edges[j].neighbour->edges.size(); k++) {

                // find edge with neighbour == our cell, then set its boundary also true
                if (cells[i].edges[j].neighbour->edges[k].neighbour == &cells[i]) {

                    // set that boundary true
                    cells[i].edges[j].neighbour->edges[k].is_boundary = true;
                }
            }
        }
    }
}


// loads structure and places it in meshpoints
// structures are given by a list of points specifying the verticies
// verticies need to be in clockwise orientation!!!
template <typename CellType>
int Mesh<CellType>::add_struct(vector<Point>* pts, double dist_a, double safety_dist, string structname) {
    
    // load structpoints from file
    vector<Point> verticies;
    ifstream file("../src/files/" + structname + ".csv");
    string line;

    while (getline(file, line)) {
        stringstream ss(line);
        string x_string, y_string;

        if (getline(ss, x_string, ',') && getline(ss, y_string, ',')) {
            double x = stod(x_string);
            double y = stod(y_string);
            verticies.emplace_back(x, y);
        }
    }
    file.close();

    // calculate corresponding seeds in inner_seeds/outer_seeds
    vector<Point> edge_vectors;
    vector<Point> a = verticies;
    vector<Point> b;
    vector<Point> midpoints;
    for (int i = 0; i < verticies.size(); i++) {
        b.emplace_back(verticies[(i+1)%verticies.size()]);
    }
    for (int i = 0; i < verticies.size(); i++) {
        edge_vectors.emplace_back((b[i].x - a[i].x), (b[i].y - a[i].y));
        midpoints.emplace_back((a[i].x + b[i].x)/2.0, (a[i].y + b[i].y)/2.0);
    }
    vector<Point> normal_vectors;
    for (int i = 0; i<edge_vectors.size(); i++) {
        normal_vectors.emplace_back((-1*(edge_vectors[i].y)/(sqrt(edge_vectors[i].x*edge_vectors[i].x + edge_vectors[i].y*edge_vectors[i].y))),
                                    ((edge_vectors[i].x)/(sqrt(edge_vectors[i].x*edge_vectors[i].x + edge_vectors[i].y*edge_vectors[i].y))));
    }
    vector<Point> inner_points;
    vector<Point> outer_points;
    for (int i = 0; i< normal_vectors.size(); i++) {
        inner_points.emplace_back((midpoints[i].x - dist_a * normal_vectors[i].x), (midpoints[i].y - dist_a * normal_vectors[i].y));
        outer_points.emplace_back((midpoints[i].x + dist_a * normal_vectors[i].x), (midpoints[i].y + dist_a * normal_vectors[i].y));
    }


    // remove seeds inside of struct (using ray-casting-algorithm, ray along x direction) and seeds to close to struct
    vector<Point> pts_removed;
    for (int i = 0; i < (*pts).size(); i++) {
        Point pt = (*pts)[i];
        int counter = 0;
        double min_dist_to_outer_seeds = 2;
        for (int j = 0; j < verticies.size(); j++) {
            double crossing_x = a[j].x + ((pt.y - a[j].y)/(b[j].y - a[j].y)) * (b[j].x - a[j].x);
            //cout << i << " _ " << (pt.y > a[j].y && pt.y < b[j].y) << " " << (pt.y < a[j].y && pt.y > b[j].y) << endl;
            if (((pt.y > a[j].y && pt.y < b[j].y) || (pt.y < a[j].y && pt.y > b[j].y)) && (pt.x < crossing_x)) {
                counter += 1;
            }

            double dist_to_outer_seed = sqrt((pt.x - outer_points[j].x)*(pt.x - outer_points[j].x) + (pt.y - outer_points[j].y)*(pt.y - outer_points[j].y));
            if (dist_to_outer_seed < min_dist_to_outer_seeds) {
                min_dist_to_outer_seeds = dist_to_outer_seed;
            }
        }

        if (counter%2 == 0 && min_dist_to_outer_seeds > safety_dist) {
            pts_removed.push_back((*pts)[i]);
            //(*pts).erase((*pts).begin() + i);
        }
    }
    *pts = pts_removed;

    // return int = inner_seeds.size(), set *pts =  {inner_seeds, outer_seeds, pts}
    inner_points.push_back(Point(0.5, 0.5));
    vector<Point> new_pts;
    for (int i = 0; i<inner_points.size(); i++) {
        new_pts.push_back(inner_points[i]);
    }
    for (int i = 0; i<outer_points.size(); i++) {
        new_pts.push_back(outer_points[i]);
    }
    for (int i = 0; i<(*pts).size(); i++) {
        new_pts.push_back((*pts)[i]);
    }
    *pts = new_pts;

    return inner_points.size();


}


// ------------------------------------------------------------------------------------------------

// step function
template <typename CellType>
double Mesh<CellType>::step_func2D(Point x, double t, Point p0, Point v, double a, double b) {
    if (x.x - v.x*t >= p0.x && x.x - v.x*t < p0.x + a && x.y - v.y*t >= p0.y && x.y - v.y*t < p0.y + b) {
        return 1;
    } else {
        return 0;
    }
}

// gaussian function
template <typename CellType>
double Mesh<CellType>::gaussian2D(Point x, double t, Point p0, Point v, double A, double sigma) {
    return A * exp(-((((x.x - v.x*t) - p0.x)*((x.x - v.x*t) - p0.x) + ((x.y - v.y*t) - p0.y)*((x.y - v.y*t) - p0.y))/(2*sigma*sigma)));
}


// legendre basis functions 2D
template <typename CellType>
double Mesh<CellType>::legendre_basisfunc2D(Point x, int n) {
    if (n == 0) {
        return 1;
    } else if (n == 1) {
        return x.x;
    } else if (n == 2) {
        return x.y;
    } else if (n == 3) {
        return x.x * x.y;
    } else if (n == 4) {
        return 0.5*(3*x.x*x.x - 1);
    } else if (n == 5) {
        return 0.5*(3*x.y*x.y - 1);
    } else
        cout << "this basisfunc does not exist (legendre_basisfunc2D)" << endl;
        return INFINITY;
}

// transform Xi to x
template <typename CellType>
Point Mesh<CellType>::ksi_to_x(Point ksi, int index) {
    double x_x = ((cells[index].edges[1].length/2)*ksi.x) + cells[index].seed.x;
    double x_y = ((cells[index].edges[0].length/2)*ksi.y) + cells[index].seed.y;
    return Point(x_x, x_y);
}

// transform x to Xi
template <typename CellType>
Point Mesh<CellType>::x_to_ksi(Point x, int index) {
    double ksi_x = (2/(cells[index].edges[1].length)) * (x.x - cells[index].seed.x);
    double ksi_y = (2/(cells[index].edges[0].length)) * (x.y - cells[index].seed.y);
    return Point(ksi_x, ksi_y);
}


// initalizes step function on a DG_Q_Cell grid using L2 projection and gaussian quadrature
template <typename CellType>
void Mesh<CellType>::DG_2D_initialize_step_function(Point p0, Point v, double a, double b) {

    // stuff for gaussian quadrature
    vector<double> gauss_quad2D_weights = {25.0/81.0, 40.0/81.0, 25.0/81.0, 40.0/81.0, 64.0/81.0, 40.0/81.0, 25.0/81.0, 40.0/81.0, 25.0/81.0};
    vector<Point> gauss_quad2D_points = {Point(-sqrt(3.0/5.0), -sqrt(3.0/5.0)), Point(0, -sqrt(3.0/5.0)), Point(sqrt(3.0/5.0), -sqrt(3.0/5.0)), Point(-sqrt(3.0/5.0), 0), Point(0, 0), Point(sqrt(3.0/5.0), 0), Point(-sqrt(3.0/5.0), sqrt(3.0/5.0)), Point(0, sqrt(3.0/5.0)), Point(sqrt(3.0/5.0), sqrt(3.0/5.0))};
    vector<double> one_over_M_jj = {1.0/4.0, 3.0/4.0, 3.0/4.0, 9.0/4.0, 5.0/4.0, 5.0/4.0};

    // go through all cells    
    for (int k = 0; k < cells.size(); k++) {
        // go through all coefficients
        for (int j = 0; j < cells[k].Q.size(); j++) {
            
            // calc correct coeff using gauss quad
            double sum = 0;
            for (int i = 0; i<gauss_quad2D_weights.size(); i++) {
                sum += one_over_M_jj[j] * step_func2D(ksi_to_x(gauss_quad2D_points[i], k), 0, p0, v, a, b) * legendre_basisfunc2D(gauss_quad2D_points[i], j) * gauss_quad2D_weights[i];
            }
            
            // set coeff
            cells[k].Q[j] = sum;
        }
    }
}


// initalizes gaussian function on a DG_Q_Cell grid using L2 projection and gaussian quadrature
template <typename CellType>
void Mesh<CellType>::DG_2D_initialize_gaussian_function(Point p0, double A, double sigma) {

    // stuff for gaussian quadrature
    vector<double> gauss_quad2D_weights = {25.0/81.0, 40.0/81.0, 25.0/81.0, 40.0/81.0, 64.0/81.0, 40.0/81.0, 25.0/81.0, 40.0/81.0, 25.0/81.0};
    vector<Point> gauss_quad2D_points = {Point(-sqrt(3.0/5.0), -sqrt(3.0/5.0)), Point(0, -sqrt(3.0/5.0)), Point(sqrt(3.0/5.0), -sqrt(3.0/5.0)), Point(-sqrt(3.0/5.0), 0), Point(0, 0), Point(sqrt(3.0/5.0), 0), Point(-sqrt(3.0/5.0), sqrt(3.0/5.0)), Point(0, sqrt(3.0/5.0)), Point(sqrt(3.0/5.0), sqrt(3.0/5.0))};
    vector<double> one_over_M_jj = {1.0/4.0, 3.0/4.0, 3.0/4.0, 9.0/4.0, 5.0/4.0, 5.0/4.0};

    // go through all cells    
    for (int k = 0; k < cells.size(); k++) {
        // go through all coefficients
        for (int j = 0; j < cells[k].Q.size(); j++) {
            
            // calc correct coeff using gauss quad
            double sum = 0;
            for (int i = 0; i<gauss_quad2D_weights.size(); i++) {
                sum += one_over_M_jj[j] * gaussian2D(ksi_to_x(gauss_quad2D_points[i], k), 0, p0, Point(0.5, 0.5), A, sigma) * legendre_basisfunc2D(gauss_quad2D_points[i], j) * gauss_quad2D_weights[i];
            }
            
            // set coeff
            cells[k].Q[j] = sum;
        }
    }
}



// SAVE DATA TO FILE ------------------------------------------------------------------------------
// saves mesh into csv file readable for python script
template <typename CellType>
void Mesh<CellType>::save_mesh(string folder_name, bool cartesian, int N_row, int sim_order, int boundary_cond, bool is_repeating, double total_sim_time, string addon, int counter, double t_sim) {

    // make sure folders exist
    string path = "src/files";
    string folder_path = path + "/" + folder_name;
    if (!filesystem::exists(path)) {filesystem::create_directories(path);}
    if (!filesystem::exists(folder_path)) {filesystem::create_directories(folder_path);}

    // open file
    string filename;
    
    string t_sim_s = to_string(int(total_sim_time));
    string t_sim_ds = to_string(int(total_sim_time * 10) % 10);

    filename = folder_path + "/" + (cartesian ? "c" : "v") + "_n" + to_string(N_row) + "_" + (is_same_v<CellType, Euler_Cell> ? "FV" : "DG") + to_string(sim_order) + "_BC" + to_string(boundary_cond) + (is_repeating ? "p" : "") + "_" + t_sim_s + "_" + t_sim_ds + "s_" + addon + "_step" + to_string(counter) + ".csv";

    if (counter == 0) {
        cout << "file storage format: " << filename << endl;
        cout << "\nsnap nr. : delta_t : t_sim, Time: [ELAPSED < ETA]" << endl << "---------------------------------------------------\n";
    }

    ofstream output_file(filename);
    
    // get maximum edge number for column correction later on
    int max_edge_nr = 0;
    for (int i = 0; i<cells.size(); i++) {
        if (cells[i].edges.size()>=max_edge_nr) {
            max_edge_nr = cells[i].edges.size();
        }
    }

    // save the mesh in the following format: ax1, ay1, bx1, by1 ; ax2, ... ; ..; | Q
    output_file << "seed.x, seed.y | a.x, a.y ; a.x ... ; | Quantities" << endl;

    for (int i = 0; i<cells.size(); i++) {

        output_file << cells[i].seed.x << "," << cells[i].seed.y << "|";

        for (int j = 0; j<cells[i].edges.size(); j++) {

            output_file << cells[i].edges[j].a.x << ','
                        << cells[i].edges[j].a.y
                        << ";";
        }

        // correct for empty columns such that further values are always at the same column
        for (int a = 0; a<(max_edge_nr-cells[i].edges.size()); a++) {
            output_file << ";";
        }

        output_file << "|";

        // store sim time
        output_file << t_sim << ",";

        // if cell type is Euler_cell save rho, u, v, E, P, gamma
        if constexpr(is_same_v<CellType, Euler_Cell>) {
            output_file << cells[i].rho << "," << cells[i].u << "," << cells[i].v << "," << cells[i].E << "," << cells[i].get_P() << "," << cells[i].gamma;
        }

        // if cell type is DG_Q_Cell save all Q out of vector
        if constexpr(is_same_v<CellType, DG_Q_Cell>) {
            for (int l = 0; l < N_basisfunc - 1; l++) {
                output_file << cells[i].Q(l) << ",";
            }
            output_file << cells[i].Q(N_basisfunc - 1);
        }

        output_file << endl;

    }

    output_file.close();

}

// load data from file
template <typename CellType>
vector<Data> Mesh<CellType>::load_file(string file_name) {
    
    ifstream file(file_name);
    vector<Data> data;
    
    if (!file.is_open()) {
        cerr << "Error opening file: " << file_name << endl;
        return data;
    }
    
    string line;
    getline(file, line); // skip header
    
    while (getline(file, line)) {
        stringstream ss(line);
        string seed_data, coords_data, q_data;
        
        if (!getline(ss, seed_data, '|') ||
            !getline(ss, coords_data, '|') ||
            !getline(ss, q_data, '|')) {
            continue;
        }
        
        Data entry;
        
        // parse seed points
        stringstream seed_ss(seed_data);
        string value;
        while (getline(seed_ss, value, ',')) {
            entry.seed_points.push_back(stod(value));
        }
        
        // parse polygons
        stringstream coords_ss(coords_data);
        string segment;
        while (getline(coords_ss, segment, ';')) {
            if (!segment.empty()) {
                vector<double> point;
                stringstream point_ss(segment);
                while (getline(point_ss, value, ',')) {
                    point.push_back(stod(value));
                }
                entry.polygon.push_back(point);
            }
        }
        if (!entry.polygon.empty()) {
            entry.polygon.push_back(entry.polygon[0]);
        }
        
        // parse quantities
        stringstream q_ss(q_data);
        while (getline(q_ss, value, ',')) {
            entry.quantities.push_back(stod(value));
        }
        
        data.push_back(entry);
    }
    
    file.close();

    return data;
}

// load mesh from file (FV + Voronoi only)
template <typename CellType>
void Mesh<CellType>::load_mesh_from_file(string file_name, bool is_repeating) {

    vector<Data> data;
    data = load_file(file_name);

    // generate mesh based on seedpoints
    vector<Point> pts;
    pts.reserve(data.size());
    for (int i = 0; i < data.size(); i++) {
        pts.push_back(Point(data[i].seed_points[0], data[i].seed_points[1]));
    }
    this->generate_vmesh2D(pts, is_repeating);

    // load quantity values into mesh
    for (int i = 0; i < cells.size(); i++) {
        cells[i].rho = data[i].quantities[1];
        cells[i].u = data[i].quantities[2];
        cells[i].v = data[i].quantities[3];
        cells[i].E = data[i].quantities[4];
        cells[i].gamma = data[i].quantities[6];
        cells[i].index = i;
    }
    is_cartesian = false;
    
    cout << "loaded IC from file: " << file_name << endl;
}


// ------------------------------------------------------------------------------------------------
// calc timestep using CFL for euler eq
template <typename CellType>
double Mesh<CellType>::dt_CFL_euler(double CFL) {
    
    double min_dt = 1;

    #pragma omp parallel for
    for (int i = 0; i < cells.size(); i++) {

        double c_i = sqrt(cells[i].get_gamma() * (cells[i].get_P()/cells[i].get_rho()));
        double R_i = sqrt(cells[i].volume / M_PI); // 2D -> hence circle
        double v_abs = sqrt((cells[i].get_u()* cells[i].get_u()) + (cells[i].get_u()* cells[i].get_u()));

        double dt_i = CFL * (R_i/(c_i + v_abs));
        if (dt_i< min_dt) {
            min_dt = dt_i;
        }
    }

    return min_dt;
}

template <typename CellType>
Mesh<CellType>::~Mesh() {}