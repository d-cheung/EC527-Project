#include "IntegralImage.h"
#include "bitmap_image.hpp"
#include <cmath>

const float IntegralImage::cR = (float).2989;
const float IntegralImage::cG = (float).5870;
const float IntegralImage::cB = (float).1140;


float IntegralImage::getValue(int y, int x)
{
	return this->Matrix[y][x];
}

void IntegralImage::setValue(int y, int x, float value)
{ 
	this->Matrix[y][x] = value;
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

	this->Matrix = new float*[height];
	for (int ii = 0; ii < height; ii++)
		Matrix[ii] = new float[width];

}

IntegralImage::~IntegralImage()
{
	if (this->Matrix != NULL) {
		for (int ii = 0; ii < Height; ii++)
			delete [] Matrix[ii];
		delete [] Matrix;
	}
	
}

//Loop Unrolling w/ Multiple Accumulators
IntegralImage * IntegralImage::FromImage(bitmap_image &image)
{
	IntegralImage * pic = new IntegralImage(image.width(), image.height());

	float rowsum0 = 0;
	float rowsum1 = 0;
	float rowsum2 = 0;
	float rowsum3 = 0;
	unsigned char red0, green0, blue0;
	unsigned char red1, green1, blue1;
	unsigned char red2, green2, blue2;
	unsigned char red3, green3, blue3;
	unsigned int x;

	for (x = 0; x < image.width()-4; x+=4)
	{
		image.get_pixel(x+0, 0, red0, green0, blue0);
		image.get_pixel(x+1, 0, red1, green1, blue1);
		image.get_pixel(x+2, 0, red2, green2, blue2);
		image.get_pixel(x+3, 0, red3, green3, blue3);
		rowsum0 += (cR * red0 + cG * green0 + cB * blue0) / (float)255;
		rowsum1 = (cR * red1 + cG * green1 + cB * blue1) / (float)255;
		rowsum2 = (cR * red2 + cG * green2 + cB * blue2) / (float)255;
		rowsum3 = (cR * red3 + cG * green3 + cB * blue3) / (float)255;
		pic->setValue(0, x+0, rowsum0);
		rowsum1 += rowsum0;
		pic->setValue(0, x+1, rowsum1);
		rowsum2 += rowsum1;
		pic->setValue(0, x+2, rowsum2);
		rowsum3 += rowsum2;
		pic->setValue(0, x+3, rowsum3);
		rowsum0 = rowsum3;
	}
	for (; x < image.width(); x++)
	{
		image.get_pixel(x+0, 0, red0, green0, blue0);
		rowsum0 += (cR * red0 + cG * green0 + cB * blue0) / (float)255;
		pic->setValue(0, x+0, rowsum0);
	}



	for (unsigned int y = 1; y < image.height(); y++)
	{
		rowsum0 = 0;
		rowsum1 = 0;
		rowsum2 = 0;
		rowsum3 = 0;
		for (x = 0; x < image.width()-4; x+=4)
		{
			image.get_pixel(x+0, y, red0, green0, blue0);
			image.get_pixel(x+1, y, red1, green1, blue1);
			image.get_pixel(x+2, y, red2, green2, blue2);
			image.get_pixel(x+3, y, red3, green3, blue3);
			rowsum0 += (cR * red0 + cG * green0 + cB * blue0) / (float)255;
			rowsum1 = (cR * red1 + cG * green1 + cB * blue1) / (float)255;
			rowsum2 = (cR * red2 + cG * green2 + cB * blue2) / (float)255;
			rowsum3 = (cR * red3 + cG * green3 + cB * blue3) / (float)255;

			// integral image is rowsum + value above		
			pic->setValue(y, x+0, rowsum0 + pic->getValue(y - 1, x+0));
			rowsum1 += rowsum0;
			pic->setValue(y, x+1, rowsum1 + pic->getValue(y - 1, x+1));
			rowsum2 += rowsum1;
			pic->setValue(y, x+2, rowsum2 + pic->getValue(y - 1, x+2));
			rowsum3 += rowsum2;
			pic->setValue(y, x+3, rowsum3 + pic->getValue(y - 1, x+3));
			rowsum0 = rowsum3;
		}
		for(; x < image.width(); x++)
		{
			image.get_pixel(x+0, y, red0, green0, blue0);
			rowsum0 += (cR * red0 + cG * green0 + cB * blue0) / (float)255;
			pic->setValue(y, x+0, rowsum0 + pic->getValue(y - 1, x+0));
		}
	}

	return pic;
}


float IntegralImage::BoxIntegral(int row, int col, int rows, int cols)
{
	// The subtraction by one for row/col is because row/col is inclusive.
	int r1 = std::min(row, Height) - 1;
	int c1 = std::min(col, Width) - 1;
	int r2 = std::min(row + rows, Height) - 1;
	int c2 = std::min(col + cols, Width) - 1;

	float A = 0, B = 0, C = 0, D = 0;

	if (r1 >= 0 && c1 >= 0) A = Matrix[r1][c1];
	if (r1 >= 0 && c2 >= 0) B = Matrix[r1][c2];
	if (r2 >= 0 && c1 >= 0) C = Matrix[r2][c1];
	if (r2 >= 0 && c2 >= 0) D = Matrix[r2][c2];

	return std::max((float)0, A - B - C + D);
}

// Get Haar Wavelet X repsonse
float IntegralImage::HaarX(int row, int column, int size)
{
	return BoxIntegral(row - size / 2, column, size, size / 2)
	- 1 * BoxIntegral(row - size / 2, column - size / 2, size, size / 2);
}

// Get Haar Wavelet Y repsonse
float IntegralImage::HaarY(int row, int column, int size)
{
	return BoxIntegral(row, column - size / 2, size / 2, size)
	- 1 * BoxIntegral(row - size / 2, column - size / 2, size / 2, size);
}

