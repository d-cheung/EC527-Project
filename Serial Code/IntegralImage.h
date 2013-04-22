#ifndef INTEGRALIMAGE_H
#define INTEGRALIMAGE_H

#include "bitmap_image.hpp"

class IntegralImage
{
private:
	static const float cR;
	static const float cG;
	static const float cB;
	float ** Matrix;

	IntegralImage(int width, int height);

public:
	int Width, Height;

	//Get the value at [y,x]
	float getValue(int y, int x);
	
	//Set the value at [y,x]
	void setValue(int y, int x, float value);

	//get the IntegralImage from Bitmap
	static IntegralImage FromImage(bitmap_image image);
	
	//Compute the BoxIntegral
	float BoxIntegral(int row, int col, int rows, int cols);

	// Get Haar Wavelet X repsonse
	float HaarX(int row, int column, int size);

	// Get Haar Wavelet Y repsonse
	float HaarY(int row, int column, int size);

	IntegralImage();
};

#endif