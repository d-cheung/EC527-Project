/*
  To compile (assuming cuda install in /usr/local/cuda ): 
  
  nvcc -I/usr/local/cuda/include -c -o cudakernel.o cudakernel.cu

  gcc -I/usr/local/cuda/include -c -o cudademo.o cudademo.c 

  gcc -o cudademo cudademo.o cudakernel.o -L/usr/local/cuda/lib -lcuda -lcudart -lm 

*/

#include <stdio.h>
#include "cuda.h"
#include "cuda_runtime_api.h"

main(int argc, char **argv){
  
  /* registers */
  int n;

  /* device ID */
  int devid;
  
  /* device count */
  int devcount;
  
  /* number of entries in arrays */
  int N = 512;

  /* pointer to host array */
  float *h_array;
  
  /* pointer to gpu device array */
  float *g_array1, *g_array2;
  
  /* find number of device in current "context" */
  cudaGetDevice(&devid);

  /* find how many devices are available */
  cudaGetDeviceCount(&devcount);

  /* allocate array on host (via CUDA) */
  cudaMallocHost((void**) &h_array, N*sizeof(float));

  /* allocate arrays on device (via CUDA) */
  cudaMalloc((void**) &g_array1, N*sizeof(float));
  cudaMalloc((void**) &g_array2, N*sizeof(float));

  /* fill up the host array */
  for(n=0;n<N;++n)
    h_array[n] = n;
  
  /* copy from host array to device array */
  cudaMemcpy(g_array1, h_array, N*sizeof(float), cudaMemcpyHostToDevice);

  /* invoke kernel on device */
  void cudakernel(int N, float *g_data, float *g_result);
  cudakernel(N, g_array1, g_array2);

  /* copy from device array to host array */
  cudaMemcpy(h_array, g_array2, N*sizeof(float), cudaMemcpyDeviceToHost);
  
  /* print out results */
  for(n=0;n<N;++n){
    printf("%f ", h_array[n]);
  }
  printf("\n");
}
