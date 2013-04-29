#include "IntegralImage.h"
#include "bitmap_image.hpp"
#include <cmath>
#include "timespec.h"

#define CPG 0.55           // Cycles per GHz -- Adjust to your computer

const float     IntegralImage::cR = (float).2989;
const float     IntegralImage::cG = (float).5870;
const float     IntegralImage::cB = (float).1140;

float IntegralImage::getValue(int y, int x)
{
#ifdef FLOAT_SINGLE
	return this->Matrix[y*Height+x];
#endif

#ifdef FLOAT_DOUBLE
	return this->Matrix[y][x];
#endif
}

void IntegralImage::setValue(int y, int x, float value)
{ 
#ifdef FLOAT_SINGLE
	this->Matrix[y*Height+x] = value;
#endif

#ifdef FLOAT_DOUBLE
	this->Matrix[y][x] = value;
#endif
}

IntegralImage::IntegralImage()
{
	this->Width = 0;
	this->Height = 0;
	this->Matrix = NULL;
}

IntegralImage::IntegralImage(int width, int height)
{
	this->Width = width;
	this->Height = height;

#ifdef FLOAT_SINGLE
	this->Matrix = (float*)malloc(sizeof(float)*width*height);
#endif

#ifdef FLOAT_DOUBLE
	this->Matrix = (float**)malloc(height*sizeof(float *));
	for (int ii = 0; ii < height; ii++)
		(this->Matrix)[ii] = (float*)malloc(width*sizeof(float));
#endif
}

IntegralImage::~IntegralImage()
{

	if (Matrix != NULL)
	{
		float ** temp[this->Height];
		cudaMemcpy(temp, this->Matrix, (this->Height)*sizeof(float *), cudaMemcpyDeviceToHost);

		for (int ii = 0; ii < this->Height; ii++)
			cudaFree(temp[ii]);

		cudaFree(this->Matrix);
	}

}

/* example kernel */
__global__ void kernelFromImageCols(unsigned char * d_image, float ** d_pic, int Width, int Height, unsigned int row_increment, unsigned int bytes_per_pixel)
{

	float colsum = (float)(0.0);
	unsigned char red, green, blue;
  	unsigned int  row_increment_ = row_increment;
	unsigned int  bytes_per_pixel_ = bytes_per_pixel;
	int width = Width;
	int height = Height;
	const float cR = (float).2989;
	const float cG = (float).5870;
	const float cB = (float).1140;

	unsigned int x = threadIdx.x + blockDim.x*blockIdx.x;
	if (x >= width) return;

	for (unsigned int y = 0; y < height; y++)
	{
		blue  = d_image[(y * row_increment_) + (x * bytes_per_pixel_ + 0)];
		green = d_image[(y * row_increment_) + (x * bytes_per_pixel_ + 1)];
		red   = d_image[(y * row_increment_) + (x * bytes_per_pixel_ + 2)];

		colsum += (cR * red + cG * green + cB * blue) / (float)255;
		d_pic[y][x] = colsum;
	}
  
}

__global__ void kernelFromImageRows(float ** d_pic, int Width, int Height){

	int width = Width;
	int height = Height;

	int y = threadIdx.x + blockDim.x*blockIdx.x;
	if (y >= height) return;

	float rowsum = d_pic[y][0];

	for (unsigned int x = 1; x < width; x++)
	{
		rowsum += d_pic[y][x];
		d_pic[y][x] = rowsum;
	}
}

IntegralImage * IntegralImage::FromImage(bitmap_image &h_image)
{
	int devid;
	int devcount;
	cudaGetDevice(&devid);
	cudaGetDeviceCount(&devcount);

	/* choose 256 threads per block for high occupancy */
	int ThreadsPerBlock = 256;
	int BlocksPerGrid;

	IntegralImage * h_pic = new IntegralImage();
	h_pic->Width  = h_image.width();
	h_pic->Height = h_image.height();

#ifdef FLOAT_DOUBLE
	cudaMalloc((void **) &(h_pic->Matrix), (h_pic->Height)*sizeof(float *));
	float ** temp[h_pic->Height];
	cudaMemcpy(temp, h_pic->Matrix, (h_pic->Height)*sizeof(float *), cudaMemcpyDeviceToHost);
	for (int ii = 0; ii < h_pic->Height; ii++)
		cudaMalloc((void **) &(temp[ii]), (h_pic->Width)*sizeof(float));
	cudaMemcpy(h_pic->Matrix, temp, (h_pic->Height)*sizeof(float *), cudaMemcpyHostToDevice);

#endif

#ifdef FLOAT_SINGLE
	cudaMallocPitch((void **) &(h_pic->Matrix), &(h_pic->Pitch), (h_pic->Width)*sizeof(float), h_pic->Height);
#endif


/*
	float * d_pic;
	cudaMalloc((void **) &d_pic, (h_pic->Width * h_pic->Height)*sizeof(float));
*/

	unsigned char * d_image;
	cudaMalloc((void**) &d_image, (h_image.length_)*sizeof(unsigned char));
	cudaMemcpy(d_image, h_image.data_, (h_image.length_)*sizeof(unsigned char), cudaMemcpyHostToDevice);

/*
	int clock_gettime(clockid_t clk_id, struct timespec *tp);

	struct timespec time1, time2;
	struct timespec time_stamp;
	clock_gettime(CLOCK_REALTIME, &time1);
*/

	/* find number of blocks */
	BlocksPerGrid = ((h_pic->Width / ThreadsPerBlock) + 1);
//	BlocksPerGrid = (((h_pic->Width < ThreadsPerBlock) || (h_pic->Width % ThreadsPerBlock)) ? ((h_pic->Width / ThreadsPerBlock) + 1) : (h_pic->Width / ThreadsPerBlock));
	kernelFromImageCols <<< BlocksPerGrid, ThreadsPerBlock >>> (d_image, h_pic->Matrix, h_pic->Width, h_pic->Height, h_image.row_increment_, h_image.bytes_per_pixel_);
 	cudaDeviceSynchronize();

	BlocksPerGrid = ((h_pic->Height / ThreadsPerBlock) + 1);
//	BlocksPerGrid = (((h_pic->Height < ThreadsPerBlock) || (h_pic->Height % ThreadsPerBlock)) ? ((h_pic->Height / ThreadsPerBlock) + 1) : (h_pic->Height / ThreadsPerBlock));
	kernelFromImageRows <<< BlocksPerGrid, ThreadsPerBlock >>> (h_pic->Matrix, h_pic->Width, h_pic->Height);
 	cudaDeviceSynchronize();

/*
	clock_gettime(CLOCK_REALTIME, &time2);
	time_stamp = diff(time1,time2);
	printf("Computation only: %ld\n", (long int)((double)(CPG)*(double)(GIG * time_stamp.tv_sec + time_stamp.tv_nsec)));
*/

	/* copy from device array to host array */
//	cudaMemcpy(h_pic->Matrix, d_pic, (h_image.width()*h_image.height())*sizeof(float), cudaMemcpyDeviceToHost);


//	cudaFree(d_pic);
	cudaFree(d_image);
	return h_pic; 
}



float IntegralImage::BoxIntegral(int row, int col, int rows, int cols)
{
	// The subtraction by one for row/col is because row/col is inclusive.
	int r1 = std::min(row, Height) - 1;
	int c1 = std::min(col, Width) - 1;
	int r2 = std::min(row + rows, Height) - 1;
	int c2 = std::min(col + cols, Width) - 1;

	float A = 0, B = 0, C = 0, D = 0;

#ifdef FLOAT_SINGLE
	if (r1 >= 0 && c1 >= 0) A = Matrix[r1*Height+c1];
	if (r1 >= 0 && c2 >= 0) B = Matrix[r1*Height+c2];
	if (r2 >= 0 && c1 >= 0) C = Matrix[r2*Height+c1];
	if (r2 >= 0 && c2 >= 0) D = Matrix[r2*Height+c2];
#endif
#ifdef FLOAT_DOUBLE
	if (r1 >= 0 && c1 >= 0) A = Matrix[r1][c1];
	if (r1 >= 0 && c2 >= 0) B = Matrix[r1][c2];
	if (r2 >= 0 && c1 >= 0) C = Matrix[r2][c1];
	if (r2 >= 0 && c2 >= 0) D = Matrix[r2][c2];
#endif

	return std::max((float)0, A - B - C + D);
}


/// <summary>
/// Get Haar Wavelet X repsonse
/// </summary>
/// <param name="row"></param>
/// <param name="column"></param>
/// <param name="size"></param>
/// <returns></returns>
float IntegralImage::HaarX(int row, int column, int size)
{
	return BoxIntegral(row - size / 2, column, size, size / 2)
	- 1 * BoxIntegral(row - size / 2, column - size / 2, size, size / 2);
}

/// <summary>
/// Get Haar Wavelet Y repsonse
/// </summary>
/// <param name="row"></param>
/// <param name="column"></param>
/// <param name="size"></param>
/// <returns></returns>
float IntegralImage::HaarY(int row, int column, int size)
{
	return BoxIntegral(row, column - size / 2, size / 2, size)
	- 1 * BoxIntegral(row - size / 2, column - size / 2, size / 2, size);
}


