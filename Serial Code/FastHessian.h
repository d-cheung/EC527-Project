#ifndef FASTHESSIAN_H
#define FASTHESSIAN_H

#include "IntegralImage.h"

class FastHessian
{

public:
	static std::vector<IPoint> * FastHessian::getIpoints(float thresh, int octaves, int init_sample, IntegralImage img);
    FastHessian(float thresh, int octaves, int init_sample, IntegralImage img);
    std::vector<IPoint> * getIpoints();

private: 
    /// Reponse Layer 
    class ResponseLayer
    {
	
		public:

		int width, height, step, filter;
		float * responses;
		char * laplacian;

		ResponseLayer()
		{
		}

		ResponseLayer(int width, int height, int step, int filter)
		{
			this->width = width;
			this->height = height;
			this->step = step;
			this->filter = filter;

			responses = new float[width * height];
			laplacian = new char[width * height];
		}

		~ResponseLayer()
		{
			delete [] responses;
			delete [] laplacian;
		}

		char getLaplacian(int row, int column)
		{
			return laplacian[row * width + column];
		}

		char getLaplacian(int row, int column, ResponseLayer src)
		{
			int scale = this->width / src.width;
			return laplacian[(scale * row) * width + (scale * column)];
		}

		float getResponse(int row, int column)
		{
			return responses[row * width + column];
		}

		float getResponse(int row, int column, ResponseLayer src)
		{
			int scale = this->width / src.width;
			return responses[(scale * row) * width + (scale * column)];
		}
	};

    /// These are passed in
    float thresh;
    int octaves;
    int init_sample;
    IntegralImage img;

    /// These get built
    std::vector<IPoint> * ipts;
    std::vector<ResponseLayer> * responseMap;

    void buildResponseMap();
    void buildResponseLayer(ResponseLayer rl);
    bool isExtremum(int r, int c, ResponseLayer t, ResponseLayer m, ResponseLayer b);
    void interpolateExtremum(int r, int c, ResponseLayer t, ResponseLayer m, ResponseLayer b);
    double * BuildDerivative(int r, int c, ResponseLayer t, ResponseLayer m, ResponseLayer b);
    double ** BuildHessian(int r, int c, ResponseLayer t, ResponseLayer m, ResponseLayer b);
	double ** Inverse(double ** m);
	double * MMM_neg_3x3_3x1(double ** A, double * B);
};

#endif