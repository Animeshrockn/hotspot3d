
#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include"3D.h"



float hotspot_stencil_core(float temp_center, float temp_top, float temp_bottom, float temp_west,float temp_east,float temp_north, float temp_south, float power_center ) {
    #pragma HLS inline off
    #ifndef _SYNTHESIS_
    // if (check == 0 && i ==0 && l ==0 && k ==0){
        
     //printf("temp_center = %f \n ", temp_center);
    //printf("temp_top  = %f  \n ", temp_top);
     //printf("temp_bottom  = %f \n    ", temp_bottom);
    //printf("temp_west = %f    \n ", temp_west);
    //printf("temp_east = %f    \n ", temp_east);
    //printf("temp_north = %f   \n  ", temp_north);
     printf("temp_south = %f   \n  ", temp_south);
    // printf("power_center = %f \n    ", power_center);
    
    #endif

    //center * cc + north * cn + south * cs + east * ce + west * cw + top * ct + bottom * cb + (dt / Cap) * pIn[z* GRID_ROWS * GRID_COLS  + r * GRID_COLS + c] + ct * AMB_TEMP;
    float tmp = (float)temp_center * cc;
	float tmp0 = (float)tmp + (float)temp_north * cn;
	float tmp1 = (float)tmp0 + (float)temp_south * cs;
    float tmp2 = (float)tmp1 + (float)temp_east * ce;
	float tmp3 = (float)tmp2 + (float)temp_west * cw;
    float tmp4 = (float)tmp3 + (float)temp_top * ct;
	float tmp5 = (float)tmp4 + (float)temp_bottom * cb;
    float tmp6 = (float)tmp5 + (float)power_center * stepDivCap;
    float tmp7 = (float)tmp6 + (float)AMB_TEMP * ct;
	return tmp7;
}

void hotspot3D(float power[TILE_LAYERS * GRID_ROWS * GRID_COLS], float temp[(TILE_LAYERS + 2) * (GRID_ROWS)*GRID_COLS], float result[TILE_LAYERS * GRID_ROWS * GRID_COLS], int l)
{
    int r, i, j, k, ii;
    float temp_top[PARA_FACTOR], temp_bottom[PARA_FACTOR], temp_west[PARA_FACTOR], temp_east[PARA_FACTOR], 
          temp_center[PARA_FACTOR], power_center[PARA_FACTOR], temp_north[PARA_FACTOR], temp_south[PARA_FACTOR];

    float temp_rf_T [PARA_FACTOR][GRID_COLS* 2 / PARA_FACTOR + 1];
    // #pragma HLS array_partition variable=temp_rf_T complete dim=0

    float temp_rf [PARA_FACTOR][GRID_COLS* 2 / PARA_FACTOR + 1];
    // #pragma HLS array_partition variable=temp_rf complete dim=0

    float temp_rf_B [PARA_FACTOR][GRID_COLS* 2 / PARA_FACTOR + 1];
    // #pragma HLS array_partition variable=temp_rf_B complete dim=0
    
    /*
        00  01  02  03  04  05  06  07  08  09        0  1  2  3  4  5  6  7  8  9  10   
        10  11  12  13  14  15  16  17  18  19        00 02 04 06 08 00 02 04 06 08 10  
        20	21	22	23	24	25	26	27	28	29        01 03 05 07 09 01 03 05 07 09 11
        30	31	32	33	34	35	36	37	38	39        c = grid_cols/para_factor = 10/2 = 5        
        40	41	42	43	44	45	46	47	48	49  ->    element in prev row = [k][0]
        50	51	52	53	54	55	56	57	58	59        element in next row = [k][GRID_COLS* 2 / PARA_FACTOR]         
        60	61	62	63	64	65	66	67	68	69                            = [k][10]
        70	71	72	73	74	75	76	77	78	79
        80	81	82	83	84	85	86	87	88	89
        90	91	92	93	94	95	96	97	98	99

    
    */

    //Loop for copying first row 
    for ( i = 0 ; i < GRID_COLS/PARA_FACTOR; i++) {
        #pragma HLS pipeline II=1
        for (ii = 0; ii < PARA_FACTOR; ii++) {
            #pragma HLS unroll
            temp_rf[ii][i] = temp[GRID_COLS*GRID_ROWS + i*PARA_FACTOR + ii];            
        }
    }

    for (  j=0; i < GRID_COLS * 2 / PARA_FACTOR + 1; i++, j++) {
        #pragma HLS pipeline II=1
        for (ii = 0; ii < PARA_FACTOR; ii++) {
            #pragma HLS unroll
            temp_rf_T[ii][i] = temp[j*PARA_FACTOR + ii];
            temp_rf[ii][i] = temp[GRID_COLS*GRID_ROWS + j*PARA_FACTOR + ii];
            temp_rf_B[ii][i] = temp[2*GRID_COLS*GRID_ROWS + j*PARA_FACTOR + ii];            
        }
    }


    for (i = 0; i < GRID_COLS / PARA_FACTOR * GRID_ROWS  ; i++) {
        #pragma HLS pipeline II=1
        for (k = 0; k < PARA_FACTOR; k++) {
            #pragma HLS unrolls
            temp_center[k] = temp_rf[k][GRID_COLS / PARA_FACTOR];

            temp_top[k] = (l==0)?temp_center[k]:temp_rf_T[k][GRID_COLS / PARA_FACTOR];
            temp_bottom[k] =  temp_rf_B[k][GRID_COLS / PARA_FACTOR];
            temp_north[k] = temp_rf[k][0];
            temp_south[k] = (i >= GRID_COLS / PARA_FACTOR * (GRID_ROWS - 1) && (i == GRID_COLS / PARA_FACTOR * GRID_ROWS -1)) ? temp_center[k] : temp_rf[k][GRID_COLS / PARA_FACTOR * 2];


            temp_west[k] = ((i % (GRID_COLS / PARA_FACTOR)) == 0 && k == 0) ? temp_center[k] : temp_rf[(k - 1 + PARA_FACTOR) % PARA_FACTOR][GRID_COLS / PARA_FACTOR - (k == 0) ];
            temp_east[k] = ((i % (GRID_COLS / PARA_FACTOR)) == (GRID_COLS / PARA_FACTOR - 1) && k == PARA_FACTOR - 1) ? temp_center[k] : temp_rf[(k + 1 + PARA_FACTOR) % PARA_FACTOR][GRID_COLS / PARA_FACTOR + (k == (PARA_FACTOR - 1)) ];
            

            power_center[k] = power[i * PARA_FACTOR + k];
            //float temp_center, float temp_top, float temp_bottom, float temp_west,float temp_east,float temp_north, float temp_south, float power_center)
            result[i * PARA_FACTOR + k] = hotspot_stencil_core(temp_center[k], temp_top[k], temp_bottom[k], temp_west[k], temp_east[k], temp_north[k], temp_south[k], power_center[k]);
        }
       if(i<5)
       {
       printf("\nTEMP_RF before shift FOR i=%d \n",i);
        for (k = 0; k < PARA_FACTOR; k++) {
            for (j = 0; j < GRID_COLS * 2 / PARA_FACTOR + 1; j++) {
                printf(" %f ",temp_rf[k][j]);
            }
             printf("\n");
        }
       }
        for (k = 0; k < PARA_FACTOR; k++) {
            #pragma hls unroll
            for (j = 0; j < GRID_COLS * 2 / PARA_FACTOR; j++) {
                #pragma hls unroll
                temp_rf_T[k][j] = temp_rf_T[k][j + 1];
                temp_rf[k][j] = temp_rf[k][j + 1];
                temp_rf_B[k][j] = temp_rf_B[k][j + 1];
            }        
            temp_rf_T[k][GRID_COLS * 2 / PARA_FACTOR] = temp[ GRID_COLS * 2 + (i+1) * PARA_FACTOR + k];
            temp_rf[k][GRID_COLS * 2 / PARA_FACTOR] = temp[GRID_COLS*GRID_ROWS + GRID_COLS * 2 + (i+1 ) * PARA_FACTOR + k - GRID_COLS];
            temp_rf_B[k][GRID_COLS * 2 / PARA_FACTOR] = temp[2*GRID_COLS*GRID_ROWS + GRID_COLS * 2 + (i+1 ) * PARA_FACTOR + k];
        }

       if(i<5)
       {
        printf("\nTEMP_RF after shift FOR i=%d \n",i);
        for (k = 0; k < PARA_FACTOR; k++) {
            for (j = 0; j < GRID_COLS * 2 / PARA_FACTOR + 1; j++) {
                printf(" %f ",temp_rf[k][j]);
            }
             printf("\n");
        }
       }
    }

}

void load(float *temp_inner,float *power_inner, float *tempIn, float *powerIn, int l)
{
  memcpy(temp_inner, tempIn + ((l-1) * GRID_ROWS * GRID_COLS), sizeof(float) * (TILE_LAYERS+2) * (GRID_ROWS)*GRID_COLS);
  memcpy(power_inner, powerIn + (l * GRID_ROWS * GRID_COLS), sizeof(float) * TILE_LAYERS * GRID_ROWS * GRID_COLS);
}

void store(float *tempOut, float *result_inner, int l)
{
  memcpy(tempOut + l * GRID_ROWS * GRID_COLS, result_inner, sizeof(float) * TILE_LAYERS * GRID_ROWS * GRID_COLS);
}

extern "C" {
    void workload(float powerIn[GRID_ROWS * GRID_COLS * GRID_LAYERS],
              float tempIn[GRID_ROWS * GRID_COLS * GRID_LAYERS],
              float tempOut[GRID_ROWS * GRID_COLS * GRID_LAYERS]
              )
{

#pragma HLS INTERFACE m_axi port = tempIn offset = slave bundle = gmem
#pragma HLS INTERFACE m_axi port = powerIn offset = slave bundle = gmem
#pragma HLS INTERFACE m_axi port = tempOut offset = slave bundle = gmem

#pragma HLS INTERFACE s_axilite port = tempIn bundle = control
#pragma HLS INTERFACE s_asizexilite port = powerIn bundle = control
#pragma HLS INTERFACE s_axilite port = tempOut bundle = control

#pragma HLS INTERFACE s_axilite port = return bundle = control

    float result_inner[TILE_LAYERS * GRID_ROWS * GRID_COLS];
    float temp_inner[(TILE_LAYERS + 2) * (GRID_ROWS)*GRID_COLS];
    float power_inner[TILE_LAYERS* GRID_ROWS * GRID_COLS];
    int i;

    ITER_LOOP:for ( i = 0; i < ITERATIONS ; i++)
        {
        int l;
        TILE_LOOP:for (l = 0; l < GRID_LAYERS; l+= TILE_LAYERS)
        {
            
            load( temp_inner,power_inner, tempIn, powerIn,  l);

            hotspot3D(power_inner, temp_inner, result_inner, l);
            
            store(tempOut, result_inner, l);
           
        }
        float *temp = tempIn;
        tempIn = tempOut;
        tempOut = temp; 
    } 
    return;
}
}