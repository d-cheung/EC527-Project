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
	return this->Matrix[y*Height+x];
//	return this->Matrix[y][x];
}

void IntegralImage::setValue(int y, int x, float value)
{ 
	this->Matrix[y*Height+x] = value;
//	this->Matrix[y][x] = value;
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

	this->Matrix = (float*)malloc(sizeof(float)*width*height);

/*
	this->Matrix = new float*[height];
	for (int ii = 0; ii < height; ii++)
		Matrix[ii] = new float[width];
*/

}

IntegralImage::~IntegralImage()
{
/*
	if (this->Matrix != NULL) {
		for (int ii = 0; ii < Height; ii++)
			delete [] Matrix[ii];
		delete [] Matrix;
	}
*/
	delete [] Matrix;

}

/* example kernel */
__global__ void kernelFromImageRows(unsigned char * d_image, float * d_pic, int Width, int Height, unsigned int row_increment, unsigned int bytes_per_pixel)
{

	float colsum = 0;
	unsigned char red, green, blue;
  	unsigned int  row_increment_ = row_increment;
	unsigned int  bytes_per_pixel_ = bytes_per_pixel;
	int width = Width;
	int height = Height;
	const float cR = (float).2989;
	const float cG = (float).5870;
	const float cB = (float).1140;

	unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
	if (y >= height) return;

	for (unsigned int x = 0; x < width; x++)
	{
		blue  = d_image[(y * row_increment_) + (x * bytes_per_pixel_ + 0)];
		green = d_image[(y * row_increment_) + (x * bytes_per_pixel_ + 1)];
		red   = d_image[(y * row_increment_) + (x * bytes_per_pixel_ + 2)];

		colsum += (cR * red + cG * green + cB * blue) / (float)255;
		d_pic[y*height + x] = colsum;
	}
  
}

__global__ void kernelFromImageCols(float * d_pic, int Width, int Height){

	int width = Width;
	int height = Height;

	int x = threadIdx.x + blockDim.x*blockIdx.x;
	if (x >= width) return;

	float rowsum = d_pic[x];

	for (int y = 1; y < height; y++)
	{
		rowsum += d_pic[y*height + x];
		d_pic[y*height + x] = rowsum;
	}
}



IntegralImage * IntegralImage::FromImage(bitmap_image &h_image)
{
	IntegralImage * h_pic = new IntegralImage(h_image.width(), h_image.height());

	float * d_pic;
//	cudaMallocPitch((void **)&(d_pic->Matrix), &pitchF,(h_image.width())*sizeof(float),h_image.height());
	cudaMalloc((void **) &d_pic, (h_image.width() * h_image.height())*sizeof(float));

	unsigned char * d_image;
	cudaMalloc((void**) &d_image, (h_image.length_)*sizeof(unsigned char));
	cudaMemcpy(d_image, h_image.data_, (h_image.length_)*sizeof(unsigned char), cudaMemcpyHostToDevice);

	int devid;
	int devcount;
	cudaGetDevice(&devid);
	cudaGetDeviceCount(&devcount);

	/* choose 256 threads per block for high occupancy */
	int ThreadsPerBlock = 256;
	int BlocksPerGrid;

	int clock_gettime(clockid_t clk_id, struct timespec *tp);

	struct timespec time1, time2;
	struct timespec time_stamp;
	clock_gettime(CLOCK_REALTIME, &time1);


	/* find number of blocks */
	BlocksPerGrid = (h_image.height() % ThreadsPerBlock) ? ((h_image.height() % ThreadsPerBlock) + 1) : (h_image.height() % ThreadsPerBlock);
	kernelFromImageRows <<< BlocksPerGrid, ThreadsPerBlock >>> (d_image, d_pic, h_image.width(), h_image.height(), h_image.row_increment_, h_image.bytes_per_pixel_);
 	cudaThreadSynchronize();

	BlocksPerGrid = (h_image.width() % ThreadsPerBlock) ? ((h_image.width() % ThreadsPerBlock) + 1) : (h_image.width() % ThreadsPerBlock);
	kernelFromImageCols <<< BlocksPerGrid, ThreadsPerBlock >>> (d_pic, h_image.width(), h_image.height());
 	cudaThreadSynchronize();

	clock_gettime(CLOCK_REALTIME, &time2);
	time_stamp = diff(time1,time2);
	printf("Computation only: %ld\n", (long int)((double)(CPG)*(double)(GIG * time_stamp.tv_sec + time_stamp.tv_nsec)));

	/* copy from device array to host array */
//	cudaMemcpy2D(h_pic->Matrix,pitchF,d_pic->Matrix, pitchF,(h_image.width())*sizeof(float),h_image.height(),cudaMemcpyDeviceToHost);
	cudaMemcpy(h_pic->Matrix, d_pic, (h_image.width()*h_image.height())*sizeof(float), cudaMemcpyDeviceToHost);


	cudaFree(d_pic);
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
	if (r1 >= 0 && c1 >= 0) A = Matrix[r1*Height+c1];
	if (r1 >= 0 && c2 >= 0) B = Matrix[r1*Height+c2];
	if (r2 >= 0 && c1 >= 0) C = Matrix[r2*Height+c1];
	if (r2 >= 0 && c2 >= 0) D = Matrix[r2*Height+c2];
/*
	if (r1 >= 0 && c1 >= 0) A = Matrix[r1][c1];
	if (r1 >= 0 && c2 >= 0) B = Matrix[r1][c2];
	if (r2 >= 0 && c1 >= 0) C = Matrix[r2][c1];
	if (r2 >= 0 && c2 >= 0) D = Matrix[r2][c2];
*/

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

