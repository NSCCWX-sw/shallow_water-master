#ifndef GRID_H
#define GRID_H
#define size_grid 258
#define MINI_GRID_WIDTH 256
#define MINI_CELL_WIDTH 32
#define MINI_GRID_NUM 8
#include <netcdfcpp.h>

#define LEFT  0
#define RIGHT 1
#define UP    2
#define DOWN  3

// testing
#define rc 1

const double gravity = 9.81;
enum bathy { flat, xslope, yslope };

extern double** h0_global;
extern double** h_global;
extern double** u_global;
extern double** v_global;
extern double** zeta_global;
extern double** zeta_star_global;

class Grid{
  private:
    double x_length; // length in x direction
    double y_length; // length in y direction
    double** h0; // the depth in 2 dimension
    double** h; // the depth in 2 dimension
    double** u; // the depth in 2 dimension
    double** v; // the depth in 2 dimension
    double** zeta; // the depth in 2 dimension
    double** zeta_star; // the depth in 2 dimension
    double hmax; // maximum depth
    double hmin; // minimum depth
    //grid parameters
    double dx; // grid spacing in x direction
    double dy; // grid spacing in y direction
    bathy slope; // the depth type
    bool initialized; // the depth type
    int grid_dimx, grid_dimy;
    
  public:
    Grid(double, double, double, double, bathy);
    void setDepth(double m_h);
    double get_x();
    double get_y();
    double get_dx();
    double get_dy();
    double** get_h();
    int get_grid_dimx();
    int get_grid_dimy();
    void show();
    void maxDepth();
    void minDepth();
    void initialize(int argc, char* argv[]);
    int solve(Grid myGrid, float dt,float time, float epsilon, int argc, char* argv[]);

};
#endif
