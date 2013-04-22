#ifndef INTEGRALIMAGE_H
#define INTEGRALIMAGE_H


public class IntegralImage
{
private:
	const float cR = .2989f;
	const float cG = .5870f;
	const float cB = .1140f;
	float[][] Matrix;

	IntegralImage(int width, int height);

public:
	int Width, Height;

	//Get the value at [y,x]
	float getValue(int y, int x);
	
	//Set the value at [y,x]
	void setValue(int y, int x);

	//get the IntegralImage from Bitmap
	static IntegralImage FromImage(Bitmap image);
	
	//Compute the BoxIntegral
	float BoxIntegral(int row, int col, int rows, int cols);

	// Get Haar Wavelet X repsonse
	float HaarX(int row, int column, int size);

	// Get Haar Wavelet Y repsonse
	float HaarY(int row, int column, int size);
}

#endif