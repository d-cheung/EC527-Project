#ifndef INTEGRALIMAGE_H
#define INTEGRALIMAGE_H

#include "bitmap_image.hpp"
#include <pthread.h>

class IntegralImage
{
private:
	static const int NUM_THREADS;
	pthread_barrier_t barr;

	static const float cR;
	static const float cG;
	static const float cB;
//	float * Matrix;
	float ** Matrix;

	IntegralImage(int width, int height);
	
	//PThread work function
	void FromImage_work(void * threadarg);

public:
	int Width, Height;

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

struct thread_data {
	int thread_id;
	IntegralImage * pic;
	bitmap_image &image;
};

#endif
