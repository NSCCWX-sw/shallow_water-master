/*************************************************************************
    > File Name: slave.c
    > Author: not wusihai
    > Mail: not wusihai18@gmail.com
    > Created Time: 2019年10月22日 星期二 09时22分27秒
 ************************************************************************/

#include "slave.h"
#include "register.h"
#include <math.h>
#include <simd.h>
#include <stdio.h>
#include <string.h>

#define ALLSYN() athread_syn(ARRAY_SCOPE, 0xffff)

#define dx 10
#define dy 10
#define epsilon 0.1

#define GRID_WIDTH 256

// definitions in solve
#define MINI_CELL_SIZE MINI_CELL_WIDTH* MINI_CELL_WIDTH
#define MINI_CELL_WIDTH 32
#define MINI_CELL_NUM 20
#define MINI_GRID_SIZE MINI_GRID_WIDTH* MINI_GRID_WIDTH
#define MINI_GRID_WIDTH 256
#define MINI_GRID_NUM 8

// testing
#define rc 1

typedef struct {
    double* p;
    int* y;
    int* x;
    int* t;
} par;

// MPI boundary
extern double** h0_global;
extern double** h_global;
extern double** u_global;
extern double** v_global;
extern double** zeta_global;
extern double** zeta_star_global;

extern double** u_lefter_global;
extern double** u_righter_global;
extern double** u_upper_global;
extern double** u_lower_global;

extern double** v_lefter_global;
extern double** v_righter_global;
extern double** v_upper_global;
extern double** v_lower_global;

extern double** zeta_star_lefter_global;
extern double** zeta_star_righter_global;
extern double** zeta_star_upper_global;
extern double** zeta_star_lower_global;

extern double** zeta_lefter_global;
extern double** zeta_righter_global;
extern double** zeta_upper_global;
extern double** zeta_lower_global;

// rc boundary
extern double** rc_u_lefter_global;
extern double** rc_u_righter_global;
extern double** rc_v_lefter_global;
extern double** rc_v_righter_global;
extern double** rc_zeta_lefter_global;
extern double** rc_zeta_righter_global;
extern double** rc_zeta_star_lefter_global;
extern double** rc_zeta_star_righter_global;
extern double** rc_h_lefter_global;
extern double** rc_h_righter_global;
extern double** rc_h0_lefter_global;
extern double** rc_h0_righter_global;


__thread_local double u_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
                    v_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
                    zeta_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
                    zeta_star_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
                    h_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
                    h0_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2];


void func_setDepth_1(par* ptr)
{
    // int i, j;
    // volatile double h0_slave[MINI_CELL_WIDTH+2][MINI_CELL_WIDTH+2];
    // volatile unsigned int get_reply, put_reply;
    // int my_id, offset_edges, edge_x, edge_y;
    // my_id = athread_get_id(-1);
    // offset_edges = my_id * CELL_SIZE;
    // for (j = 0; j < CELL_SIZE; j++) {
    //     edge_y = (offset_edges + j) / GRID_WIDTH;
    //     edge_x = (offset_edges + j) % GRID_WIDTH;
    //     if (edge_y == 0 || edge_y == 257 || edge_x == 0 && edge_x == 257) {
    //         h0_slave[j] = 0;
    //     } else {
    //         h0_slave[j] = *(ptr->p);
    //     }
    // }
    // put_reply = 0;
    // athread_put(PE_MODE, &h0_slave[0],
    //     &h0_global[offset_edges / GRID_WIDTH][offset_edges % GRID_WIDTH],
    //     CELL_SIZE * 8, &put_reply, 0, 0);
    // while (put_reply != 1)
    //     ;
}

// mode 1 initialization
void func_initialize(par* ptr)
{
    // int edge_x, edge_y, offset_edges;
    // int i, j;
    // double arr_slave[CELL_SIZE];
    // volatile int my_id = athread_get_id(-1);
    // volatile unsigned int get_reply, put_reply;
    // get_reply = 0;
    // int peak = HALF_GRID_WIDTH * GRID_WIDTH + HALF_GRID_WIDTH;
    // int total_size = CELL_SIZE * 8;
    // volatile int position_master = my_id * CELL_SIZE;
    // athread_get(
    //     PE_MODE,
    //     &h0_global[position_master / GRID_WIDTH][position_master % GRID_WIDTH],
    //     &arr_slave[0], total_size, &get_reply, 0, 0, 0);
    // while (get_reply != 1)
    //     ;
    // put_reply = 0;
    // // h_global
    // athread_put(
    //     PE_MODE, &arr_slave[0],
    //     &h_global[position_master / GRID_WIDTH][position_master % GRID_WIDTH],
    //     total_size, &put_reply, 0, 0);
    // while (put_reply != 1)
    //     ;
    // for (j = 0; j < CELL_SIZE; j++) {
    //     arr_slave[j] = 0;
    // }
    // put_reply = 0;
    // // zeta_star_global
    // athread_put(PE_MODE, &arr_slave[0],
    //     &zeta_star_global[position_master / GRID_WIDTH]
    //                      [position_master % GRID_WIDTH],
    //     total_size, &put_reply, 0, 0);
    // // u_global
    // athread_put(
    //     PE_MODE, &arr_slave[0],
    //     &u_global[position_master / GRID_WIDTH][position_master % GRID_WIDTH],
    //     total_size, &put_reply, 0, 0);
    // // v_global
    // athread_put(
    //     PE_MODE, &arr_slave[0],
    //     &v_global[position_master / GRID_WIDTH][position_master % GRID_WIDTH],
    //     total_size, &put_reply, 0, 0);
    // while (put_reply != 3)
    //     ;

    // // zeta_global
    // put_reply = 0;
    // for (j = 0; j < CELL_SIZE; j++) {
    //     if (((my_id * CELL_SIZE + j) / GRID_WIDTH == HALF_GRID_WIDTH || (my_id * CELL_SIZE + j) / GRID_WIDTH == HALF_GRID_WIDTH + 1) && ((my_id * CELL_SIZE + j) % GRID_WIDTH == HALF_GRID_WIDTH || (my_id * CELL_SIZE + j) % GRID_WIDTH == HALF_GRID_WIDTH + 1)) {
    //         arr_slave[j] = 1.0;
    //     }
    // }
    // athread_put(
    //     PE_MODE, &arr_slave[0],
    //     &zeta_global[position_master / GRID_WIDTH][position_master % GRID_WIDTH],
    //     total_size, &put_reply, 0, 0);
    // while (put_reply != 1)
    //     ;
}

void func_solve(par* ptr)
{
    int i, j, t;
    // volatile double u_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
    //     v_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
    //     zeta_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
    //     zeta_star_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
    //     h_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
    //     h0_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2];
    for(i=0; i<MINI_CELL_WIDTH + 2; i++){
        for(j=0; j<MINI_CELL_WIDTH + 2; j++){
            u_slave[j][i] = 0;
            v_slave[j][i] = 0;
            zeta_slave[j][i] = 0;
            zeta_star_slave[j][i] = 0;
            h_slave[j][i] = 0;
            h0_slave[j][i] = 0;
        }
    }
    volatile doublev4 boundary_in, boundary_out;
    volatile double convert[MINI_CELL_WIDTH], temp_left[MINI_CELL_WIDTH], temp_right[MINI_CELL_WIDTH];
    volatile double he, hw, hn, hs;
    volatile int position_x, position_y;
    volatile int my_id = athread_get_id(-1);
    volatile int row_id, col_id;
    GET_ROW(row_id);
    GET_COL(col_id);
    volatile unsigned int get_reply, put_reply;
    int total_size = (MINI_CELL_WIDTH + 2) * 8;
    float dt = 0.1;
    double gravity = 9.81;
    volatile int targetID[8] = { 7, 0, 1, 2, 3, 4, 5, 6 };
    volatile int targetID_1[8] = { 1, 2, 3, 4, 5, 6, 7, 0 };
    
	position_y = my_id / MINI_GRID_NUM * MINI_CELL_WIDTH;
    position_x = my_id % MINI_GRID_NUM * MINI_CELL_WIDTH;

    double rc_u_lefter[MINI_CELL_WIDTH];
    double rc_u_righter[MINI_CELL_WIDTH];
    double rc_v_lefter[MINI_CELL_WIDTH];
    double rc_v_righter[MINI_CELL_WIDTH];
    
    t = *(ptr->t);
#ifdef rc 
        if (t == 0) {
            get_reply = 0;
            for (j = 0; j < MINI_CELL_WIDTH + 2; j++) {
                athread_get(PE_MODE, &u_global[position_y + j][position_x], &u_slave[j][0], total_size, &get_reply, 0, 0, 0);
                athread_get(PE_MODE, &zeta_global[position_y + j][position_x], &zeta_slave[j][0], total_size, &get_reply, 0, 0, 0);
            }
            while (get_reply != 2 * (MINI_CELL_WIDTH + 2))
                ;
        }else{
            get_reply = 0;
            athread_get(PE_MODE, &u_global[position_y][position_x], &u_slave[0][0], total_size, &get_reply, 0, 0, 0);
            athread_get(PE_MODE, &u_global[position_y+MINI_CELL_WIDTH+1][position_x], &u_slave[MINI_CELL_WIDTH+1][0], total_size, &get_reply, 0, 0, 0);
            athread_get(PE_MODE, &v_global[position_y][position_x], &v_slave[0][0], total_size, &get_reply, 0, 0, 0);
            athread_get(PE_MODE, &v_global[position_y+MINI_CELL_WIDTH+1][position_x], &v_slave[MINI_CELL_WIDTH+1][0], total_size, &get_reply, 0, 0, 0);
            total_size = MINI_CELL_WIDTH * 8;
            athread_get(PE_MODE, &rc_u_lefter_global[col_id][row_id*MINI_CELL_WIDTH], &rc_u_lefter[0], total_size, &get_reply, 0, 0, 0);
            athread_get(PE_MODE, &rc_u_righter_global[col_id][row_id*MINI_CELL_WIDTH], &rc_u_righter[0], total_size, &get_reply, 0, 0, 0);
            athread_get(PE_MODE, &rc_v_lefter_global[col_id][row_id*MINI_CELL_WIDTH], &rc_v_lefter[0], total_size, &get_reply, 0, 0, 0);
            athread_get(PE_MODE, &rc_v_righter_global[col_id][row_id*MINI_CELL_WIDTH], &rc_v_righter[0], total_size, &get_reply, 0, 0, 0);
            while (get_reply != 8)
                ;
            for (j = 1; j < MINI_CELL_WIDTH + 1; j++) {
                u_slave[j-1][1] = rc_u_lefter[j-1];
                v_slave[j-1][1] = rc_v_lefter[j-1];
                u_slave[j-1][MINI_CELL_WIDTH] = rc_u_righter[j-1];
                v_slave[j-1][MINI_CELL_WIDTH] = rc_v_righter[j-1];
            }
 
        }
        
#endif

#ifdef rc 
        // targetID[0] = 7;
        // targetID[1] = 0;
        // targetID[2] = 1;
        // targetID[3] = 2;
        // targetID[4] = 3;
        // targetID[5] = 4;
        // targetID[6] = 5;
        // targetID[7] = 6;
        // zeta_global[j][k+1]
        for (j = 0; j < MINI_GRID_NUM; j++) {
            boundary_out = simd_set_doublev4(zeta_slave[j * 4 + 1][1], zeta_slave[j * 4 + 2][1],zeta_slave[j * 4 + 3][1], zeta_slave[j * 4 + 4][1]);
            REG_PUTR(boundary_out, targetID[col_id]);
            REG_GETR(boundary_in);
            simd_store(boundary_in, &convert[j * 4]);
        }
        if (col_id != 7) {
            for (j = 0; j < MINI_GRID_NUM; j++) {
                zeta_slave[j * 4 + 1][MINI_CELL_WIDTH + 1] = convert[j * 4];
                zeta_slave[j * 4 + 2][MINI_CELL_WIDTH + 1] = convert[j * 4 + 1];
                zeta_slave[j * 4 + 3][MINI_CELL_WIDTH + 1] = convert[j * 4 + 2];
                zeta_slave[j * 4 + 4][MINI_CELL_WIDTH + 1] = convert[j * 4 + 3];
            }
        }
#endif
        
        get_reply = 0;
        for (j = 0; j < MINI_CELL_WIDTH + 2; j++) {
            athread_get(PE_MODE, &v_global[position_y + j][position_x], &v_slave[j][0], total_size, &get_reply, 0, 0, 0);
        }
        while (get_reply != MINI_CELL_WIDTH + 2)
                ;

        // calculation for u
        for (j = 1; j < MINI_CELL_WIDTH + 1; j++) {
            for (i = 1; i < MINI_CELL_WIDTH + 1; i++) {
                u_slave[j][i] = u_slave[j][i] - dt * gravity * (zeta_slave[j][i + 1] - zeta_slave[j][i]) / dx;
            }
        }


            
#ifdef rc 
            // zeta_slave for v_slave
            // targetID[0] = 7;
            // targetID[1] = 0;
            // targetID[2] = 1;
            // targetID[3] = 2;
            // targetID[4] = 3;
            // targetID[5] = 4;
            // targetID[6] = 5;
            // targetID[7] = 6;
            // zeta_global[j+1][k]
            for (j = 0; j < MINI_GRID_NUM; j++) {
                boundary_out = simd_set_doublev4(zeta_slave[1][j * 4 + 1], zeta_slave[1][j * 4 + 2],
                    zeta_slave[1][j * 4 + 3], zeta_slave[1][j * 4 + 4]);
                REG_PUTC(boundary_out, targetID[row_id]);
                REG_GETC(boundary_in);
                simd_store(boundary_in, &convert[j * 4]);
            }
            if (row_id != 7) {
                for (j = 0; j < MINI_GRID_NUM; j++) {
                    zeta_slave[MINI_CELL_WIDTH + 1][j * 4 + 1] = convert[j * 4];
                    zeta_slave[MINI_CELL_WIDTH + 1][j * 4 + 2] = convert[j * 4 + 1];
                    zeta_slave[MINI_CELL_WIDTH + 1][j * 4 + 3] = convert[j * 4 + 2];
                    zeta_slave[MINI_CELL_WIDTH + 1][j * 4 + 4] = convert[j * 4 + 3];
                }
            }
#endif

        // calculation for v
        for (j = 1; j < MINI_CELL_WIDTH + 1; j++) {
            for (i = 1; i < MINI_CELL_WIDTH + 1; i++) {
                v_slave[j][i] = v_slave[j][i] - dt * gravity * (zeta_slave[j + 1][i] - zeta_slave[j][i]) / dy;
            }
        }

#ifdef rc 
            // u_slave
            // targetID[0] = 1;
            // targetID[1] = 2;
            // targetID[2] = 3;
            // targetID[3] = 4;
            // targetID[4] = 5;
            // targetID[5] = 6;
            // targetID[6] = 7;
            // targetID[7] = 0;
            // u[j][k-1]
            for (j = 0; j < MINI_GRID_NUM; j++) {
                boundary_out = simd_set_doublev4(u_slave[j * 4 + 1][MINI_CELL_WIDTH],
                    u_slave[j * 4 + 2][MINI_CELL_WIDTH],
                    u_slave[j * 4 + 3][MINI_CELL_WIDTH],
                    u_slave[j * 4 + 4][MINI_CELL_WIDTH]);
                REG_PUTR(boundary_out, targetID_1[col_id]);
                REG_GETR(boundary_in);
                simd_store(boundary_in, &convert[j * 4]);
            }
            if (col_id != 0) {
                for (j = 0; j < MINI_GRID_NUM; j++) {
                    u_slave[j * 4 + 1][0] = convert[j * 4];
                    u_slave[j * 4 + 2][0] = convert[j * 4 + 1];
                    u_slave[j * 4 + 3][0] = convert[j * 4 + 2];
                    u_slave[j * 4 + 4][0] = convert[j * 4 + 3];
                }
            }

            // v_slave[j-1][k]
            // targetID[0] = 1;
            // targetID[1] = 2;
            // targetID[2] = 3;
            // targetID[3] = 4;
            // targetID[4] = 5;
            // targetID[5] = 6;
            // targetID[6] = 7;
            // targetID[7] = 0;
            for (j = 0; j < MINI_GRID_NUM; j++) {
                boundary_out = simd_set_doublev4(v_slave[MINI_CELL_WIDTH][4 * j + 1],
                    v_slave[MINI_CELL_WIDTH][4 * j + 2],
                    v_slave[MINI_CELL_WIDTH][4 * j + 3],
                    v_slave[MINI_CELL_WIDTH][4 * j + 4]);
                REG_PUTC(boundary_out, targetID_1[row_id]);
                REG_GETC(boundary_in);
                simd_store(boundary_in, &convert[j * 4]);
           }
            if (row_id != 0) {
                for (j = 0; j < MINI_GRID_NUM; j++) {
                    v_slave[0][4 * j + 1] = convert[j * 4];
                    v_slave[0][4 * j + 2] = convert[j * 4 + 1];
                    v_slave[0][4 * j + 3] = convert[j * 4 + 2];
                    v_slave[0][4 * j + 4] = convert[j * 4 + 3];
                }
             }

            // zeta_star[j][k-1]
            // targetID[0] = 1;
            // targetID[1] = 2;
            // targetID[2] = 3;
            // targetID[3] = 4;
            // targetID[4] = 5;
            // targetID[5] = 6;
            // targetID[6] = 7;
            // targetID[7] = 0;
            for (j = 0; j < MINI_GRID_NUM; j++) {
                boundary_out = simd_set_doublev4(zeta_star_slave[4 * j + 1][MINI_CELL_WIDTH],
                    zeta_star_slave[4 * j + 2][MINI_CELL_WIDTH],
                    zeta_star_slave[4 * j + 3][MINI_CELL_WIDTH],
                    zeta_star_slave[4 * j + 4][MINI_CELL_WIDTH]);
                REG_PUTR(boundary_out, targetID_1[col_id]);
                REG_GETR(boundary_in);
                simd_store(boundary_in, &convert[4 * j]);
            }
            if (col_id != 0) {
                for (j = 0; j < MINI_GRID_NUM; j++) {
                    zeta_star_slave[j * 4 + 1][0] = convert[4 * j];
                    zeta_star_slave[j * 4 + 2][0] = convert[4 * j + 1];
                    zeta_star_slave[j * 4 + 3][0] = convert[4 * j + 2];
                    zeta_star_slave[j * 4 + 4][0] = convert[4 * j + 3];
                }
            }

            // zeta_star[j][k+1]
            // targetID[0] = 7;
            // targetID[1] = 0;
            // targetID[2] = 1;
            // targetID[3] = 2;
            // targetID[4] = 3;
            // targetID[5] = 4;
            // targetID[6] = 5;
            // targetID[7] = 6;
            for (j = 0; j < MINI_GRID_NUM; j++) {
                boundary_out = simd_set_doublev4(zeta_star_slave[4 * j + 1][1],
                    zeta_star_slave[4 * j + 2][1],
                    zeta_star_slave[4 * j + 3][1],
                    zeta_star_slave[4 * j + 4][1]);
                REG_PUTR(boundary_out, targetID[col_id]);
                REG_GETR(boundary_in);
                simd_store(boundary_in, &convert[j * 4]);
            }
            if (col_id != 7) {
                for (j = 0; j < MINI_GRID_NUM; j++) {
                    zeta_star_slave[j * 4 + 1][MINI_CELL_WIDTH + 1] = convert[j * 4];
                    zeta_star_slave[j * 4 + 2][MINI_CELL_WIDTH + 1] = convert[j * 4 + 1];
                    zeta_star_slave[j * 4 + 3][MINI_CELL_WIDTH + 1] = convert[j * 4 + 2];
                    zeta_star_slave[j * 4 + 4][MINI_CELL_WIDTH + 1] = convert[j * 4 + 3];
                }
            }      
        
            // zeta_star[j-1][k]
            // targetID[0] = 1;
            // targetID[1] = 2;
            // targetID[2] = 3;
            // targetID[3] = 4;
            // targetID[4] = 5;
            // targetID[5] = 6;
            // targetID[6] = 7;
            // targetID[7] = 0;
            for (j = 0; j < MINI_GRID_NUM; j++) {
                boundary_out = simd_set_doublev4(zeta_star_slave[MINI_CELL_WIDTH][4 * j + 1],
                    zeta_star_slave[MINI_CELL_WIDTH][4 * j + 2],
                    zeta_star_slave[MINI_CELL_WIDTH][4 * j + 3],
                    zeta_star_slave[MINI_CELL_WIDTH][4 * j + 4]);
                REG_PUTC(boundary_out, targetID_1[row_id]);
                REG_GETC(boundary_in);
                simd_store(boundary_in, &convert[j * 4]);
            }
            if (row_id != 0) {
                for (j = 0; j < MINI_GRID_NUM; j++) {
                    zeta_star_slave[0][j * 4 + 1] = convert[j * 4];
                    zeta_star_slave[0][j * 4 + 2] = convert[j * 4 + 1];
                    zeta_star_slave[0][j * 4 + 3] = convert[j * 4 + 2];
                    zeta_star_slave[0][j * 4 + 4] = convert[j * 4 + 3];
                }
            }

            // zeta_star[j+1][k]
            // targetID[0] = 7;
            // targetID[1] = 0;
            // targetID[2] = 1;
            // targetID[3] = 2;
            // targetID[4] = 3;
            // targetID[5] = 4;
            // targetID[6] = 5;
            // targetID[7] = 6;
            for (j = 0; j < MINI_GRID_NUM; j++) {
                boundary_out = simd_set_doublev4(
                    zeta_star_slave[1][4 * j + 1], zeta_star_slave[1][4 * j + 2],
                    zeta_star_slave[1][4 * j + 3], zeta_star_slave[1][4 * j + 4]);
                REG_PUTC(boundary_out, targetID[row_id]);
                REG_GETC(boundary_in);
                simd_store(boundary_in, &convert[4 * j]);
            }
            if (row_id != 7) {
                for (j = 0; j < MINI_GRID_NUM; j++) {
                    zeta_star_slave[MINI_CELL_WIDTH + 1][4 * j + 1] = convert[4 * j];
                    zeta_star_slave[MINI_CELL_WIDTH + 1][4 * j + 2] = convert[4 * j + 1];
                    zeta_star_slave[MINI_CELL_WIDTH + 1][4 * j + 3] = convert[4 * j + 2];
                    zeta_star_slave[MINI_CELL_WIDTH + 1][4 * j + 4] = convert[4 * j + 3];
                }
            }
#endif

   
#ifdef rc 
            if(t == 1999){
                // put the data back to the master processor
                total_size = MINI_CELL_WIDTH * 8;
                put_reply = 0;
                for (j = 1; j < MINI_CELL_WIDTH + 1; j++) {
                    athread_put(PE_MODE, &u_slave[j][1], &u_global[position_y + j][position_x + 1], total_size, &put_reply, 0, 0);
                    athread_put(PE_MODE, &v_slave[j][1], &v_global[position_y + j][position_x + 1], total_size, &put_reply, 0, 0);
                }
                while (put_reply != 2 * MINI_CELL_WIDTH)
                        ;    
            }else{
                // u_slave and v_slave
                // upper and lower
                put_reply = 0;
                athread_put(PE_MODE, &u_slave[1][1], &u_global[position_y + 1][position_x + 1], total_size, &put_reply, 0, 0);
                athread_put(PE_MODE, &v_slave[1][1], &v_global[position_y + 1][position_x + 1], total_size, &put_reply, 0, 0);
                athread_put(PE_MODE, &u_slave[MINI_CELL_WIDTH][1], &u_global[position_y + MINI_CELL_WIDTH][position_x + 1], total_size, &put_reply, 0, 0);
                athread_put(PE_MODE, &v_slave[MINI_CELL_WIDTH][1], &v_global[position_y + MINI_CELL_WIDTH][position_x + 1], total_size, &put_reply, 0, 0);
                while (put_reply != 4)
                        ;
                // lefer and righter
                for (j = 1; j < MINI_CELL_WIDTH + 1; j++) {
                    rc_u_lefter[j-1] = u_slave[j][1];
                    rc_v_lefter[j-1] = v_slave[j][1];
                    rc_u_righter[j-1] = u_slave[j][MINI_CELL_WIDTH];
                    rc_v_righter[j-1] = v_slave[j][MINI_CELL_WIDTH];
                }

                put_reply = 0;
                athread_put(PE_MODE, rc_u_righter, &rc_u_righter_global[my_id % MINI_GRID_NUM][position_y], total_size, &put_reply, 0, 0);
                athread_put(PE_MODE, rc_u_lefter, &rc_u_lefter_global[col_id][row_id*MINI_CELL_WIDTH], total_size, &put_reply, 0, 0);
                athread_put(PE_MODE, rc_v_righter, &rc_v_righter_global[my_id % MINI_GRID_NUM][position_y], total_size, &put_reply, 0, 0);
                athread_put(PE_MODE, rc_v_lefter, &rc_v_lefter_global[my_id % MINI_GRID_NUM][position_y], total_size, &put_reply, 0, 0);
                while (put_reply != 4)
                        ;
            }
            
#else
            // put the data back to the master processor
            total_size = MINI_CELL_WIDTH * 8;
            put_reply = 0;
            for (j = 1; j < MINI_CELL_WIDTH + 1; j++) {
                athread_put(PE_MODE, &u_slave[j][1], &u_global[position_y + j][position_x + 1], total_size, &put_reply, 0, 0);
                athread_put(PE_MODE, &v_slave[j][1], &v_global[position_y + j][position_x + 1], total_size, &put_reply, 0, 0);
            }
            while (put_reply != 2 * MINI_CELL_WIDTH)
                    ;          
#endif
}

void func_solve_1(par* ptr)
{
    int i, j, t;
    // volatile double u_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
    //     v_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
    //     zeta_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
    //     zeta_star_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
    //     h_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2],
    //     h0_slave[MINI_CELL_WIDTH + 2][MINI_CELL_WIDTH + 2];
    volatile doublev4 boundary_in, boundary_out;
    volatile double convert[MINI_CELL_WIDTH], temp_left[MINI_CELL_WIDTH], temp_right[MINI_CELL_WIDTH];
    volatile double he, hw, hn, hs;
    volatile int position_x, position_y;
    volatile int my_id = athread_get_id(-1);
    volatile int row_id, col_id;
    GET_ROW(row_id);
    GET_COL(col_id);
    volatile unsigned int get_reply, put_reply;
    int total_size = (MINI_CELL_WIDTH + 2) * 8;
    float dt = 0.1;
    double gravity = 9.81;
    volatile int targetID[8] = { 7, 0, 1, 2, 3, 4, 5, 6 };
    volatile int targetID_1[8] = { 1, 2, 3, 4, 5, 6, 7, 0 };
    
	position_y = my_id / MINI_GRID_NUM * MINI_CELL_WIDTH;
    position_x = my_id % MINI_GRID_NUM * MINI_CELL_WIDTH;

    get_reply = 0;
    for (j = 0; j < MINI_CELL_WIDTH + 2; j++) {
        athread_get(PE_MODE, &zeta_star_global[position_y + j][position_x], &zeta_star_slave[j][0], total_size, &get_reply, 0, 0, 0);
        athread_get(PE_MODE, &h_global[position_y + j][position_x], &h_slave[j][0], total_size, &get_reply, 0, 0, 0);
        athread_get(PE_MODE, &h_global[position_y + j][position_x], &h0_slave[j][0], total_size, &get_reply, 0, 0, 0);
    }
    while (get_reply != 3 * (MINI_CELL_WIDTH + 2))
        ;
    // set value for he, hw, hn, hs
    // calculate for zeta_star_slave and zeta_slave
    for (j = 1; j < MINI_CELL_WIDTH + 1; j++) {
        for (i = 1; i < MINI_CELL_WIDTH + 1; i++) {
            if (u_slave[j][i] > 0) {
                he = h_slave[j][i];
            } else if (u_slave[j][i] < 0) {
                he = h_slave[j][i + 1];
            }
            if (u_slave[j][i - 1] > 0) {
                hw = h_slave[j][i - 1];
            } else if (u_slave[j][i - 1] < 0) {
                hw = h_slave[j][i];
            }
            if (v_slave[j][i] > 0) {
                hn = h_slave[j][i];
            } else if (v_slave[j][i] < 0) {
                hn = h_slave[j + 1][i];
            }
            if (v_slave[j - 1][i] > 0) {
                hs = h_slave[j - 1][i];
            } else if (v_slave[j - 1][i] < 0) {
                hs = h_slave[j][i];
            }
            zeta_star_slave[j][i] = zeta_slave[j][i] - dt * (u_slave[j][i] * he - u_slave[j][i - 1] * hw) / dx - dt * (v_slave[j][i] * hn - v_slave[j - 1][i] * hs) / dy;
        }
    }
    for (j = 1; j < MINI_CELL_WIDTH + 1; j++) {
        for (i = 1; i < MINI_CELL_WIDTH + 1; i++) {
            zeta_slave[j][i] = (1 - epsilon) * zeta_star_slave[j][i] + 0.25 * epsilon * (zeta_star_slave[j][i - 1] + zeta_star_slave[j][i + 1] + zeta_star_slave[j - 1][i] + zeta_star_slave[j + 1][i]);
        }
    }
    for (j = 1; j < MINI_CELL_WIDTH + 1; j++) {
        for (i = 1; i < MINI_CELL_WIDTH + 1; i++) {
            h_slave[j][i] = h0_slave[j][i] + zeta_slave[j][i];
        }
    }
    total_size = MINI_CELL_WIDTH * 8;
    put_reply = 0;
    for (j = 1; j < MINI_CELL_WIDTH + 1; j++) {
        athread_put(PE_MODE, &zeta_slave[j][1], &zeta_global[position_y + j][position_x + 1], total_size, &put_reply, 0, 0);
        athread_put(PE_MODE, &h_slave[j][1], &h_global[position_y + j][position_x + 1], total_size, &put_reply, 0, 0);
        athread_put(PE_MODE, &zeta_star_slave[j][1], &zeta_star_global[position_y + j][position_x + 1], total_size, &put_reply, 0, 0);
    }
    while (put_reply != 3 * MINI_CELL_WIDTH)
            ;
}
