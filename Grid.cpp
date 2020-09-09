#include "Grid.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "mpi.h"
#include "gptl.h"

// MPI management
#define peakA_x 384
#define peakA_y 384
#define peakB_x 385
#define peakB_y 384
#define peakC_x 384
#define peakC_y 385
#define peakD_x 385
#define peakD_y 385


int my_id, numprocs;
int ndims = 2;
int dims[2] = {3, 3}; 
int periods[2] = {0, 0}; 
int coords[2], nbrs[4];
int reorder = 0;
MPI_Comm comm2d;

typedef struct {
    double* p;
    int* t;
} par;

typedef std::pair<int,int> pair;

extern "C" {
    #include <athread.h>
    void slave_func_setDepth_1(par* p);
    void slave_func_initialize();
    void slave_func_solve(par* p);
    void slave_func_solve_1(par* p);
}

int ret=GPTLinitialize();

// templte function for 2d array creation
template <class T>
T** array2d(int m, int n)
{
    T** A = new T*[n];
    T* B = new T[m * n];
    for (int i = 0; i < n; i++) {
        A[i] = &(B[i * m]);
    }
    return A;
}

double** h0_global = array2d<double>(size_grid, size_grid);
double** h_global = array2d<double>(size_grid, size_grid);
double** u_global = array2d<double>(size_grid, size_grid);
double** v_global = array2d<double>(size_grid, size_grid);
double** zeta_global = array2d<double>(size_grid, size_grid);
double** zeta_star_global = array2d<double>(size_grid, size_grid);

// MPI boundary exchange
double** u_lefter_global = array2d<double>(1, MINI_GRID_WIDTH);
double** u_righter_global = array2d<double>(1, MINI_GRID_WIDTH);
double** u_upper_global = array2d<double>(1, MINI_GRID_WIDTH);
double** u_lower_global = array2d<double>(1, MINI_GRID_WIDTH);

double** v_lefter_global = array2d<double>(1, MINI_GRID_WIDTH);
double** v_righter_global = array2d<double>(1, MINI_GRID_WIDTH);
double** v_upper_global = array2d<double>(1, MINI_GRID_WIDTH);
double** v_lower_global = array2d<double>(1, MINI_GRID_WIDTH);

double** zeta_star_lefter_global = array2d<double>(1, MINI_GRID_WIDTH);
double** zeta_star_righter_global = array2d<double>(1, MINI_GRID_WIDTH);
double** zeta_star_upper_global = array2d<double>(1, MINI_GRID_WIDTH);
double** zeta_star_lower_global = array2d<double>(1, MINI_GRID_WIDTH);

double** zeta_lefter_global = array2d<double>(1, MINI_GRID_WIDTH);
double** zeta_righter_global = array2d<double>(1, MINI_GRID_WIDTH);
double** zeta_upper_global = array2d<double>(1, MINI_GRID_WIDTH);
double** zeta_lower_global = array2d<double>(1, MINI_GRID_WIDTH);

double** h_lefter_global = array2d<double>(1, MINI_GRID_WIDTH);
double** h_righter_global = array2d<double>(1, MINI_GRID_WIDTH);
double** h_upper_global = array2d<double>(1, MINI_GRID_WIDTH);
double** h_lower_global = array2d<double>(1, MINI_GRID_WIDTH);

//rc boundary
double** rc_u_lefter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_u_righter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_v_lefter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_v_righter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_zeta_lefter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_zeta_righter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_zeta_star_lefter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_zeta_star_righter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_h_lefter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_h_righter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_h0_lefter_global = array2d<double>(MINI_GRID_WIDTH, 8);
double** rc_h0_righter_global = array2d<double>(MINI_GRID_WIDTH, 8);


// save function for octave data
// void save(double* data, int timestep,int dimx, int dimy)
void save(int timestep, int dimx, int dimy, int my_id)
{
    std::ofstream myfile;
    char filename[18];
    timestep = 1000;
    sprintf(filename, "./eta/eta%05d_%d.dat", timestep, my_id);
    myfile.open(filename);
    for (int j = 1; j < dimy-1; j++) {
        for (int i = 1; i < dimx-1; i++) {
            myfile << zeta_global[j][i] << " ";
            // myfile << h0_global[j][i] << " ";
        }
        myfile << "\n";
    }
    std::cout << "File saved successfully." << std::endl;
    myfile.close();
}

Grid::Grid(double m_x, double m_y, double m_dx, double m_dy, bathy m_slope)
    : x_length(m_x)
    , y_length(m_y)
    , dx(m_dx)
    , dy(m_dy)
    , slope(m_slope)
{
}

void Grid::setDepth(double m_h)
{
    // grid_dimx = x_length/dx+2;
    // grid_dimy = y_length/dy+2;
    /*
   * HUGE ATTENTION:
   * dimx and dimy must be in accordance with the defined size in main.cpp
   */

    grid_dimx = x_length;
    grid_dimy = y_length;
    par ptr;
    ptr.p = &m_h;

    std::cout << "dimx,dimy " << grid_dimx << "," << grid_dimy << std::endl;
    std::cout << size_grid << std::endl;

    //initialize h0
    if (slope == flat) { // flat
        std::cout << "Flat case" << std::endl;
        std::cout << grid_dimx << ", " << grid_dimy << std::endl;
        for(int j = 1; j < grid_dimy-1; j++){
          for(int i = 1; i < grid_dimx-1; i++){
            h0_global[j][i] = m_h;
          }
        }
        
        // ret = GPTLstart("assignment");
        // athread_init();
        // __real_athread_spawn((void*)slave_func_setDepth_1, &ptr);
        // athread_join();
        // ret = GPTLstop("assignment");
    } else if (slope == xslope) { // x slope
        std::cout << "X-slope" << std::endl;
        float x_step, cum_x_step;
        x_step = m_h / (grid_dimx - 3);
        cum_x_step = 0;
        for (int i = 1; i < grid_dimx - 1; i++) {
            for (int j = 1; j < grid_dimy - 1; j++) {
                h0_global[j][i] = cum_x_step;
            }
            cum_x_step += x_step;
        }
    } else if (slope == yslope) { // y slope
        std::cout << "Y-slope" << std::endl;
        float y_step, cum_y_step;
        y_step = m_h / (grid_dimy - 3);
        cum_y_step = 0;
        for (int i = 1; i < grid_dimy - 1; i++) {
            for (int j = 1; j < grid_dimx - 1; j++) {
                h0_global[i][j] = cum_y_step;
            }
            cum_y_step += y_step;
        }
    }
    // initialize edges with 0
    for (int j = 0; j < 1; j++) {
        for (int i = 0; i < grid_dimx; i++) {
            h0_global[j][i] = 0;
            h0_global[grid_dimy - j - 1][i] = 0;
        }
    }
    for (int j = 0; j < grid_dimy; j++) {
        for (int i = 0; i < 1; i++) {
            h0_global[j][i] = 0;
            h0_global[j][grid_dimx - i - 1] = 0;
        }
    }
}

double Grid::get_x() { return x_length; }
double Grid::get_y() { return y_length; }
double Grid::get_dx() { return dx; }
double Grid::get_dy() { return dy; }
double** Grid::get_h() { return h; }
int Grid::get_grid_dimx() { return grid_dimx; }
int Grid::get_grid_dimy() { return grid_dimy; }
void Grid::show()
{
    std::cout << "Grid resolution : " << x_length / dx + 2 << "x" << y_length / dy + 2 << std::endl;
    std::cout << "Max Depth : " << hmax << std::endl;
    std::cout << "Min Depth : " << hmin << std::endl;
}
void Grid::maxDepth()
{
    hmax = 0;
    for (int j = 4; j < grid_dimy - 4; j++) {
        for (int k = 4; k < grid_dimx - 4; k++) {
            hmax = fmax(hmax, h0_global[j][k]);
        }
    }
}
void Grid::minDepth()
{
    hmin = h0_global[4][4];
    for (int j = 4; j < grid_dimy - 4; j++) {
        for (int k = 4; k < grid_dimx - 4; k++) {
            hmin = fmin(hmin, h0_global[j][k]);
        }
    }
}
void Grid::initialize(int argc, char* argv[])
{
    int mode = 1;
    int time = 1000;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm2d);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Cart_coords(comm2d, my_id, 2, coords);
    MPI_Cart_shift(comm2d, 0, 1, &nbrs[2], &nbrs[3]);
    MPI_Cart_shift(comm2d, 1, 1, &nbrs[0], &nbrs[1]);
    // std::cout << my_id << ", coords" << coords[0] << " " << coords[1] << ", neighbors{l, r, u, d}" << nbrs[LEFT] << " "
    // << nbrs[RIGHT] << " " << nbrs[UP] << " " << nbrs[DOWN] << std::endl;

    if (mode == 1) {
        for(int k=0; k<grid_dimy;k++){
         for(int j=0; j<grid_dimx;j++){
           zeta_global[k][j] = 0.0;
           zeta_star_global[k][j] = 0.0;
           h_global[k][j] = h0_global[k][j]+zeta_global[k][j];
           u_global[k][j] = 0.0;
           v_global[k][j] = 0.0;
        }
        }
        // set the initial peaks    
        if (std::ceil((double)peakA_y/MINI_GRID_WIDTH)-1 == coords[0] && std::ceil((double)peakA_x/MINI_GRID_WIDTH)-1 == coords[1] ){
            zeta_global[(peakA_y-1)%MINI_GRID_WIDTH+1][(peakA_x-1)%MINI_GRID_WIDTH+1] = 1.0;
            // std::cout << "A " << my_id << " ";
            // std::cout << (peakA_y-1)%MINI_GRID_WIDTH+1 << ", " << (peakA_x-1)%MINI_GRID_WIDTH+1 << std::endl;
        }
        if (std::ceil((double)peakB_y/MINI_GRID_WIDTH)-1 == coords[0] && std::ceil((double)peakB_x/MINI_GRID_WIDTH)-1 == coords[1]){
            zeta_global[(peakB_y-1)%MINI_GRID_WIDTH+1][(peakB_x-1)%MINI_GRID_WIDTH+1] = 1.0;
            // std::cout << "B " << my_id << " ";
            // std::cout << (peakB_y-1)%MINI_GRID_WIDTH+1 << ", " << (peakB_x-1)%MINI_GRID_WIDTH+1 << std::endl;
        }
        if (std::ceil((double)peakC_x/MINI_GRID_WIDTH)-1== coords[1] && std::ceil((double)peakC_y/MINI_GRID_WIDTH)-1 == coords[0]){
            zeta_global[(peakC_y-1)%MINI_GRID_WIDTH+1][(peakC_x-1)%MINI_GRID_WIDTH+1] = 1.0;
            // std::cout << "C " << my_id << " ";
            // std::cout << (peakC_y-1)%MINI_GRID_WIDTH+1 << ", " << (peakC_x-1)%MINI_GRID_WIDTH+1 << std::endl;
        }
        if (std::ceil((double)peakD_y/MINI_GRID_WIDTH)-1 == coords[0] && std::ceil((double)peakD_x/MINI_GRID_WIDTH)-1 == coords[1]){
            zeta_global[(peakD_y-1)%MINI_GRID_WIDTH+1][(peakD_x-1)%MINI_GRID_WIDTH+1] = 1.0;
            // std::cout << "D " << my_id << " ";
            // std::cout << (peakD_y-1)%MINI_GRID_WIDTH+1 << ", " << (peakD_x-1)%MINI_GRID_WIDTH+1 << std::endl;
        }
        // __real_athread_spawn((void*)slave_func_initialize, 0);
        // athread_join();
        int count = 0;
        for (int i = 0; i < grid_dimx; i++) {
            for (int j = 0; j < grid_dimy; j++) {
                if (h_global[i][j] == 10) {
                    count += 1;
                }
            }
        }
        // std::cout << "This is the count in h_global: " << count << std::endl;
        count = 0;
        for (int i = 0; i < grid_dimx; i++) {
            for (int j = 0; j < grid_dimy; j++) {
                if (zeta_global[i][j] == 1.0) {
                    count += 1;
                }
            }
        }
        // std::cout << "This is the count in zeta_global: " << count << std::endl;
    }
    this->maxDepth();
    this->minDepth();
    initialized = 1;
    std::cout << "function initialize is done." << std::endl;
}

int Grid::solve(Grid myGrid, float dt, float time, float epsilon, int argc, char* argv[])
{   
    MPI_Request reqs[2];

    if(my_id == 0){
        std::cout << "SWW solve starts." << std::endl;
    }
    MPI_Barrier(comm2d);

    std::cout << my_id << ", coords" << coords[0] << " " << coords[1] << ", neighbors{l, r, u, d}" << nbrs[LEFT] << " "
    << nbrs[RIGHT] << " " << nbrs[UP] << " " << nbrs[DOWN] << std::endl;

    double master_bound_out[256];
    MPI_Barrier(comm2d);
    if(my_id == 0){
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "|░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░|" << std::endl;
        std::cout << "|░░░░Shallow Water 2D░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░|" << std::endl;
        std::cout << "|░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░|" << std::endl;
        std::cout << "|░░░░berkay pirlepeli░░░░2019░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░|" << std::endl;
        std::cout << "|░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░|" << std::endl;
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
    if (!initialized) {
        std::cout << "parameters are not initialized" << std::endl;
        return -1;
    } else {
        // Test stability criteria
        double min_spacing = myGrid.get_dx();
        if (min_spacing > myGrid.get_dy()) {
            min_spacing = myGrid.get_dy();
        }
        if (dt > (min_spacing / (sqrt(2 * gravity * hmax)))) {
            std::cerr << "Error! Stability Criteria is not satisfied." << std::endl;
        }
        double x = myGrid.get_x();
        double y = myGrid.get_y();
        //open netcdf file
        NcFile ofile("output.nc", NcFile::Replace);
        // add dimensions
        NcDim* dim_x = ofile.add_dim("x", grid_dimx);
        NcDim* dim_y = ofile.add_dim("y", grid_dimy);
        NcDim* dim_t = ofile.add_dim("time");

        // create variables
        NcVar* var_z = ofile.add_var("zeta", ncFloat, dim_t, dim_y, dim_x);
        var_z->add_att("units", "meter");
        var_z->add_att("long_name", "free surface height");

        double he, hw, hn, hs;
        // simulation loop
        std::cout << "Simulation Loop" << std::endl;
        par ptr;
        int time = 2000;
        // int time = 1;
        ptr.t = &time;

        int tag = 0;
        // initialize the boundaries
        for(int i=0; i<MINI_GRID_WIDTH; i++){
            u_lefter_global[0][i] = 0;
            u_righter_global[0][i] = 0;
            u_upper_global[0][i] = 0;
            u_lower_global[0][i] = 0;
            v_lefter_global[0][i] = 0;
            v_righter_global[0][i] = 0;
            v_upper_global[0][i] = 0;
            v_lower_global[0][i] = 0;
            zeta_star_lefter_global[0][i] = 0;
            zeta_star_righter_global[0][i] = 0;
            zeta_star_upper_global[0][i] = 0;
            zeta_star_lower_global[0][i] = 0;
            zeta_lefter_global[0][i] = 0;
            zeta_righter_global[0][i] = 0;
            zeta_upper_global[0][i] = 0;
            zeta_lower_global[0][i] = 0;
        }

        MPI_Status status;
        athread_init();

        ret = GPTLstart("finite difference calculation");
        for(int t = 0; t < time; t++){
            // zeta_star
            // left bound
            tag = 8;
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                master_bound_out[i] = zeta_star_global[i+1][1];
            }
            MPI_Send(master_bound_out, MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[LEFT], tag, MPI_COMM_WORLD);
            MPI_Recv(&zeta_star_righter_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[RIGHT], tag, MPI_COMM_WORLD, &status);
            if(nbrs[RIGHT] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    zeta_star_global[i+1][size_grid-1] = zeta_star_righter_global[0][i];
                }
            }
            // right bound
            tag = 9;
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                master_bound_out[i] = zeta_star_global[i+1][size_grid-2];
            }
            MPI_Send(master_bound_out, MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[RIGHT], tag, MPI_COMM_WORLD);
            MPI_Recv(&zeta_star_lefter_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[LEFT], tag, MPI_COMM_WORLD, &status);
            if(nbrs[LEFT] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    zeta_star_global[i+1][0] = zeta_star_lefter_global[0][i];
                }
            }
            // up bound
            tag = 10;
            MPI_Send(&zeta_star_global[1][1], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[UP], tag, MPI_COMM_WORLD);
            MPI_Recv(&zeta_star_lower_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[DOWN], tag, MPI_COMM_WORLD, &status);
            if(nbrs[DOWN] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    zeta_star_global[size_grid-1][i+1] = zeta_star_lower_global[0][i];
               }
            }
            // low bound
            tag = 11;
            MPI_Send(&zeta_star_global[size_grid-2][1], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[DOWN], tag, MPI_COMM_WORLD);
            MPI_Recv(&zeta_star_upper_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[UP], tag, MPI_COMM_WORLD, &status);
            if(nbrs[UP] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    zeta_star_global[0][i+1] = zeta_star_upper_global[0][i];
                }
            }

            // zeta
            // left bound
            tag = 12;
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                master_bound_out[i] = zeta_global[i+1][1];
            }
            MPI_Send(master_bound_out, MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[LEFT], tag, MPI_COMM_WORLD);
            MPI_Recv(&zeta_righter_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[RIGHT], tag, MPI_COMM_WORLD, &status);
            if(nbrs[RIGHT] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    zeta_global[i+1][size_grid-1] = zeta_righter_global[0][i];
                }
            }
            // right bound
            tag = 13;
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                master_bound_out[i] = zeta_global[i+1][size_grid-2];
            }
            MPI_Send(master_bound_out, MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[RIGHT], tag, MPI_COMM_WORLD);
            MPI_Recv(&zeta_lefter_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[LEFT], tag, MPI_COMM_WORLD, &status);
            if(nbrs[LEFT] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    zeta_global[i+1][0] = zeta_lefter_global[0][i];
                }
            }
            // up bound
            tag = 14;
            MPI_Send(&zeta_global[1][1], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[UP], tag, MPI_COMM_WORLD);
            MPI_Recv(&zeta_lower_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[DOWN], tag, MPI_COMM_WORLD, &status);
            if(nbrs[DOWN] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    zeta_global[size_grid-1][i+1] = zeta_lower_global[0][i];
                }
            }
            // low bound
            tag = 15;
            MPI_Send(&zeta_global[size_grid-2][1], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[DOWN], tag, MPI_COMM_WORLD);
            MPI_Recv(&zeta_upper_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[UP], tag, MPI_COMM_WORLD, &status);
            if(nbrs[UP] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    zeta_global[0][i+1] = zeta_upper_global[0][i];
                }
            }

#ifdef rc 
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                rc_u_lefter_global[0][i] = u_global[i+1][MINI_CELL_WIDTH*0+1];
                rc_u_lefter_global[1][i] = u_global[i+1][MINI_CELL_WIDTH*1+1];
                rc_u_lefter_global[2][i] = u_global[i+1][MINI_CELL_WIDTH*2+1];
                rc_u_lefter_global[3][i] = u_global[i+1][MINI_CELL_WIDTH*3+1];
                rc_u_lefter_global[4][i] = u_global[i+1][MINI_CELL_WIDTH*4+1];
                rc_u_lefter_global[5][i] = u_global[i+1][MINI_CELL_WIDTH*5+1];
                rc_u_lefter_global[6][i] = u_global[i+1][MINI_CELL_WIDTH*6+1];
                rc_u_lefter_global[7][i] = u_global[i+1][MINI_CELL_WIDTH*7+1];

                rc_u_righter_global[0][i] = u_global[i+1][MINI_CELL_WIDTH*1];
                rc_u_righter_global[1][i] = u_global[i+1][MINI_CELL_WIDTH*2];
                rc_u_righter_global[2][i] = u_global[i+1][MINI_CELL_WIDTH*3];
                rc_u_righter_global[3][i] = u_global[i+1][MINI_CELL_WIDTH*4];
                rc_u_righter_global[4][i] = u_global[i+1][MINI_CELL_WIDTH*5];
                rc_u_righter_global[5][i] = u_global[i+1][MINI_CELL_WIDTH*6];
                rc_u_righter_global[6][i] = u_global[i+1][MINI_CELL_WIDTH*7];
                rc_u_righter_global[7][i] = u_global[i+1][MINI_CELL_WIDTH*8];

                rc_v_lefter_global[0][i] = v_global[i+1][MINI_CELL_WIDTH*0+1];
                rc_v_lefter_global[1][i] = v_global[i+1][MINI_CELL_WIDTH*1+1];
                rc_v_lefter_global[2][i] = v_global[i+1][MINI_CELL_WIDTH*2+1];
                rc_v_lefter_global[3][i] = v_global[i+1][MINI_CELL_WIDTH*3+1];
                rc_v_lefter_global[4][i] = v_global[i+1][MINI_CELL_WIDTH*4+1];
                rc_v_lefter_global[5][i] = v_global[i+1][MINI_CELL_WIDTH*5+1];
                rc_v_lefter_global[6][i] = v_global[i+1][MINI_CELL_WIDTH*6+1];
                rc_v_lefter_global[7][i] = v_global[i+1][MINI_CELL_WIDTH*7+1];

                rc_v_righter_global[0][i] = v_global[i+1][MINI_CELL_WIDTH*1];
                rc_v_righter_global[1][i] = v_global[i+1][MINI_CELL_WIDTH*2];
                rc_v_righter_global[2][i] = v_global[i+1][MINI_CELL_WIDTH*3];
                rc_v_righter_global[3][i] = v_global[i+1][MINI_CELL_WIDTH*4];
                rc_v_righter_global[4][i] = v_global[i+1][MINI_CELL_WIDTH*5];
                rc_v_righter_global[5][i] = v_global[i+1][MINI_CELL_WIDTH*6];
                rc_v_righter_global[6][i] = v_global[i+1][MINI_CELL_WIDTH*7];
                rc_v_righter_global[7][i] = v_global[i+1][MINI_CELL_WIDTH*8];
            }
#endif

            __real_athread_spawn((void*)slave_func_solve, &ptr);
            athread_join();

#ifdef rc 
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                u_global[i+1][MINI_CELL_WIDTH*0+1] = rc_u_lefter_global[0][i];
                u_global[i+1][MINI_CELL_WIDTH*1+1] = rc_u_lefter_global[1][i];
                u_global[i+1][MINI_CELL_WIDTH*2+1] = rc_u_lefter_global[2][i];
                u_global[i+1][MINI_CELL_WIDTH*3+1] = rc_u_lefter_global[3][i];
                u_global[i+1][MINI_CELL_WIDTH*4+1] = rc_u_lefter_global[4][i];
                u_global[i+1][MINI_CELL_WIDTH*5+1] = rc_u_lefter_global[5][i];
                u_global[i+1][MINI_CELL_WIDTH*6+1] = rc_u_lefter_global[6][i];
                u_global[i+1][MINI_CELL_WIDTH*7+1] = rc_u_lefter_global[7][i];

                u_global[i+1][MINI_CELL_WIDTH*1] = rc_u_righter_global[0][i];
                u_global[i+1][MINI_CELL_WIDTH*2] = rc_u_righter_global[1][i];
                u_global[i+1][MINI_CELL_WIDTH*3] = rc_u_righter_global[2][i];
                u_global[i+1][MINI_CELL_WIDTH*4] = rc_u_righter_global[3][i];
                u_global[i+1][MINI_CELL_WIDTH*5] = rc_u_righter_global[4][i];
                u_global[i+1][MINI_CELL_WIDTH*6] = rc_u_righter_global[5][i];
                u_global[i+1][MINI_CELL_WIDTH*7] = rc_u_righter_global[6][i];
                u_global[i+1][MINI_CELL_WIDTH*8] = rc_u_righter_global[7][i];

                v_global[i+1][MINI_CELL_WIDTH*0+1] = rc_v_lefter_global[0][i];
                v_global[i+1][MINI_CELL_WIDTH*1+1] = rc_v_lefter_global[1][i];
                v_global[i+1][MINI_CELL_WIDTH*2+1] = rc_v_lefter_global[2][i];
                v_global[i+1][MINI_CELL_WIDTH*3+1] = rc_v_lefter_global[3][i];
                v_global[i+1][MINI_CELL_WIDTH*4+1] = rc_v_lefter_global[4][i];
                v_global[i+1][MINI_CELL_WIDTH*5+1] = rc_v_lefter_global[5][i];
                v_global[i+1][MINI_CELL_WIDTH*6+1] = rc_v_lefter_global[6][i];
                v_global[i+1][MINI_CELL_WIDTH*7+1] = rc_v_lefter_global[7][i];

                v_global[i+1][MINI_CELL_WIDTH*1] = rc_v_righter_global[0][i];
                v_global[i+1][MINI_CELL_WIDTH*2] = rc_v_righter_global[1][i];
                v_global[i+1][MINI_CELL_WIDTH*3] = rc_v_righter_global[2][i];
                v_global[i+1][MINI_CELL_WIDTH*4] = rc_v_righter_global[3][i];
                v_global[i+1][MINI_CELL_WIDTH*5] = rc_v_righter_global[4][i];
                v_global[i+1][MINI_CELL_WIDTH*6] = rc_v_righter_global[5][i];
                v_global[i+1][MINI_CELL_WIDTH*7] = rc_v_righter_global[6][i];
                v_global[i+1][MINI_CELL_WIDTH*8] = rc_v_righter_global[7][i];
            }
#endif

            // u
            // left bound
            tag = 0;
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                master_bound_out[i] = u_global[i+1][1];
            }
            MPI_Send(master_bound_out, MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[LEFT], tag, MPI_COMM_WORLD);
            MPI_Recv(&u_righter_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[RIGHT], tag, MPI_COMM_WORLD, &status);
            if(nbrs[RIGHT] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    u_global[i+1][size_grid-1] = u_righter_global[0][i];
                }
            }
            // right bound
            tag = 1;
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                master_bound_out[i] = u_global[i+1][size_grid-2];
            }
            MPI_Send(master_bound_out, MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[RIGHT], tag, MPI_COMM_WORLD);
            MPI_Recv(&u_lefter_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[LEFT], tag, MPI_COMM_WORLD, &status);
            if(nbrs[LEFT] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    u_global[i+1][0] = u_lefter_global[0][i];
                }
            }
            // up bound
            tag = 2;
            MPI_Send(&u_global[1][1], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[UP], tag, MPI_COMM_WORLD);
            MPI_Recv(&u_lower_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[DOWN], tag, MPI_COMM_WORLD, &status);
            if(nbrs[DOWN] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    u_global[size_grid-1][i+1] = u_lower_global[0][i];
                }
            }
            // low bound
            tag = 3;
            MPI_Send(&u_global[size_grid-2][1], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[DOWN], tag, MPI_COMM_WORLD);
            MPI_Recv(&u_upper_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[UP], tag, MPI_COMM_WORLD, &status);
            if(nbrs[UP] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    u_global[0][i+1] = u_upper_global[0][i];
                }
            }
           
            // v
            // left bound
            tag = 4;
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                master_bound_out[i] = v_global[i+1][1];
            }
            MPI_Send(master_bound_out, MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[LEFT], tag, MPI_COMM_WORLD);
            MPI_Recv(&v_righter_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[RIGHT], tag, MPI_COMM_WORLD, &status);
            if(nbrs[RIGHT] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    v_global[i+1][size_grid-1] = v_righter_global[0][i];
                }            
            }
            // right bound
            tag = 5;
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                master_bound_out[i] = v_global[i+1][size_grid-2];
            }
            MPI_Send(master_bound_out, MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[RIGHT], tag, MPI_COMM_WORLD);
            MPI_Recv(&v_lefter_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[LEFT], tag, MPI_COMM_WORLD, &status);
            if(nbrs[LEFT] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    v_global[i+1][0] = v_lefter_global[0][i];
                }
            }
            // up bound
            tag = 6;
            MPI_Send(&v_global[1][1], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[UP], tag, MPI_COMM_WORLD);
            MPI_Recv(&v_lower_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[DOWN], tag, MPI_COMM_WORLD, &status);
            if(nbrs[DOWN] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    v_global[size_grid-1][i+1] = v_lower_global[0][i];
                }
            }
            // low bound
            tag = 7;
            MPI_Send(&v_global[size_grid-2][1], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[DOWN], tag, MPI_COMM_WORLD);
            MPI_Recv(&v_upper_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[UP], tag, MPI_COMM_WORLD, &status);
            if(nbrs[UP] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    v_global[0][i+1] = v_upper_global[0][i];
                }
            }

            // h
            // left bound
            tag = 16;
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                master_bound_out[i] = h_global[i+1][1];
            }
            MPI_Send(master_bound_out, MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[LEFT], tag, MPI_COMM_WORLD);
            MPI_Recv(&h_righter_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[RIGHT], tag, MPI_COMM_WORLD, &status);
            if(nbrs[RIGHT] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    h_global[i+1][size_grid-1] = h_righter_global[0][i];
                }            
            }
            // right bound
            tag = 17;
            for(int i=0; i<MINI_GRID_WIDTH; i++){
                master_bound_out[i] = h_global[i+1][size_grid-2];
            }
            MPI_Send(master_bound_out, MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[RIGHT], tag, MPI_COMM_WORLD);
            MPI_Recv(&h_lefter_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[LEFT], tag, MPI_COMM_WORLD, &status);
            if(nbrs[LEFT] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    h_global[i+1][0] = h_lefter_global[0][i];
                }
            }
            // up bound
            tag = 18;
            MPI_Send(&h_global[1][1], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[UP], tag, MPI_COMM_WORLD);
            MPI_Recv(&h_lower_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[DOWN], tag, MPI_COMM_WORLD, &status);
            if(nbrs[DOWN] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    h_global[size_grid-1][i+1] = h_lower_global[0][i];
                }
            }
            // low bound
            tag = 19;
            MPI_Send(&h_global[size_grid-2][1], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[DOWN], tag, MPI_COMM_WORLD);
            MPI_Recv(&h_upper_global[0][0], MINI_GRID_WIDTH, MPI_DOUBLE, nbrs[UP], tag, MPI_COMM_WORLD, &status);
            if(nbrs[UP] != -1){
                for(int i=0; i<MINI_GRID_WIDTH; i++){
                    h_global[0][i+1] = h_upper_global[0][i];
                }
            }

            __real_athread_spawn((void*)slave_func_solve_1, &ptr);
            athread_join();
        }

        athread_halt();
        ret = GPTLstop("finite difference calculation");

        save(time, grid_dimx, grid_dimy, my_id);
    }
    GPTLpr(0);
    MPI_Finalize();
    return 0;
}
