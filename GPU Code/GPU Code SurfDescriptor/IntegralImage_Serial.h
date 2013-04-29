#ifndef INTEGRALIMAGE_SERIAL_H
#define INTEGRALIMAGE_SERIAL_H

#include "bitmap_image.hpp"

class IntegralImage_Serial
{
private:
	static const float cR;
	static const float cG;
	static const float cB;

	IntegralImage_Serial(int width, int height);

public:
//	float * Matrix;
	float ** Matrix;
	int Width, Height;

	//Get the value at [y,x]
	float getValue(int y, int x);
	
	//Set the value at [y,x]
	void setValue(int y, int x, float value);

	//get the IntegralImage_Serial from Bitmap
	static IntegralImage_Serial * FromImage(bitmap_image &image);
	
	//Compute the BoxIntegral
	float BoxIntegral(int row, int col, int rows, int cols);

	// Get Haar Wavelet X repsonse
	float HaarX(int row, int column, int size);

	// Get Haar Wavelet Y repsonse
	float HaarY(int row, int column, int size);

	IntegralImage_Serial();
	~IntegralImage_Serial();
};

#endif
