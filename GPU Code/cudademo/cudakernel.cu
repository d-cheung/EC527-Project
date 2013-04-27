#include "cuda.h"

/* example kernel */
__global__ void kernel(int N, float *g_data, float *g_result){

  int n = threadIdx.x + blockDim.x*blockIdx.x;
  
  /* compute squares of entries in data array */
  if(n<N)
    g_result[n] = g_data[n]*g_data[n];
  
}

/* only use extern if calling code is C */
extern "C" 
{
/* driver for kernel */

void cudakernel(int N, float *g_data, float *g_result){
  
  /* choose 256 threads per block for high occupancy */
  int ThreadsPerBlock = 256;
  
  /* find number of blocks */
  int BlocksPerGrid = (N+ThreadsPerBlock-1)/ThreadsPerBlock;
  
  /* invoke device on this block/thread grid */
  kernel <<< BlocksPerGrid, ThreadsPerBlock >>> (N, g_data, g_result);

}

}
