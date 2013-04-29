#ifndef INTEGRALIMAGE_H
#define INTEGRALIMAGE_H


#define define_min(a,b) ((a)<(b)?(a):(b))
#define define_max(a,b) ((a)>(b)?(a):(b))

#include "bitmap_image.hpp"

//#define FLOAT_SINGLE
#define FLOAT_DOUBLE


#ifdef FLOAT_SINGLE
typedef float* currentTry;
#endif

#ifdef FLOAT_DOUBLE
typedef float** currentTry;
#endif

class IntegralImage
{
private:
	IntegralImage(int width, int height);

public:
	static const float cR;
	static const float cG;
	static const float cB;
	int Width, Height;

	currentTry Matrix;
	size_t Pitch;

	//Get the value at [y,x]
	float getValue(int y, int x);

	//Set the value at [y,x]
	void setValue(int y, int x, float value);

	//get the IntegralImage from Bitmap
	static IntegralImage * FromImage(bitmap_image &image);
	
	//Compute the BoxIntegral
	float BoxIntegral(int row, int col, int rows, int cols);

	// Get Haar Wavelet X repsonse
	float HaarX(int row, int column, int size);

	// Get Haar Wavelet Y repsonse
	float HaarY(int row, int column, int size);

	IntegralImage();

	~IntegralImage();
};


__device__ inline float d_BoxIntegral(float ** img, int row, int col, int rows, int cols, int Height, int Width)
{
	int r1 = define_min(row, Height) - 1;
	int c1 = define_min(col, Width) - 1;
	int r2 = define_min(row + rows, Height) - 1;
	int c2 = define_min(col + cols, Width) - 1;

	float A = 0, B = 0, C = 0, D = 0;
	if (r1 >= 0 && c1 >= 0) A = img[r1][c1];
	if (r1 >= 0 && c2 >= 0) B = img[r1][c2];
	if (r2 >= 0 && c1 >= 0) C = img[r2][c1];
	if (r2 >= 0 && c2 >= 0) D = img[r2][c2];

	return define_max((float)0, A - B - C + D);
}

__device__ inline float d_HaarX(float ** img, int row, int column, int size, int Height, int Width)
{
	return d_BoxIntegral(img, row - size / 2,            column, size, size / 2, Height, Width)
	 - 1 * d_BoxIntegral(img, row - size / 2, column - size / 2, size, size / 2, Height, Width);
}

__device__ inline float d_HaarY(float ** img, int row, int column, int size, int Height, int Width)
{
	return d_BoxIntegral(img,            row, column - size / 2, size / 2, size, Height, Width)
	 - 1 * d_BoxIntegral(img, row - size / 2, column - size / 2, size / 2, size, Height, Width);
}


__device__ inline float d_BoxIntegral(float * img, int row, int col, int rows, int cols, size_t Pitch, int Height, int Width)
{
	int r1 = define_min(row, Height) - 1;
	int c1 = define_min(col, Width) - 1;
	int r2 = define_min(row + rows, Height) - 1;
	int c2 = define_min(col + cols, Width) - 1;
	size_t pitch = Pitch;

	float A = 0, B = 0, C = 0, D = 0;
	if (r1 >= 0 && c1 >= 0) A = *((float *)((char*)img + r1 * pitch) + c1);
	if (r1 >= 0 && c2 >= 0) B = *((float *)((char*)img + r1 * pitch) + c2);
	if (r2 >= 0 && c1 >= 0) C = *((float *)((char*)img + r2 * pitch) + c1);
	if (r2 >= 0 && c2 >= 0) D = *((float *)((char*)img + r2 * pitch) + c2);

	return define_max((float)0, A - B - C + D);
}


__device__ inline float d_HaarX(float * img, int row, int column, int size, size_t Pitch, int Height, int Width)
{
	return d_BoxIntegral(img, row - size / 2,            column, size, size / 2, Pitch, Height, Width)
	 - 1 * d_BoxIntegral(img, row - size / 2, column - size / 2, size, size / 2, Pitch, Height, Width);
}

__device__ inline float d_HaarY(float * img, int row, int column, int size, size_t Pitch, int Height, int Width)
{
	return d_BoxIntegral(img,            row, column - size / 2, size / 2, size, Pitch, Height, Width)
	 - 1 * d_BoxIntegral(img, row - size / 2, column - size / 2, size / 2, size, Pitch, Height, Width);
}



#endif
