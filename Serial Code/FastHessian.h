#ifndef FASTHESSIAN_H
#define FASTHESSIAN_H

#include "IntegralImage.h"
#include "IPoint.h"
#include <vector>

class FastHessian
{

public:
	static std::vector<IPoint> * getIpoints(float thresh, int octaves, int init_sample, IntegralImage * img);
	FastHessian(float thresh, int octaves, int init_sample, IntegralImage * img);
	~FastHessian();
	std::vector<IPoint> * getIpoints();

private: 
	// Response Layer 
	class ResponseLayer
	{
	
		public:

		int width, height, step, filter;
		float * responses;
		unsigned char * laplacian;

		ResponseLayer();
		ResponseLayer(int width, int height, int step, int filter);
		~ResponseLayer();
		unsigned char getLaplacian(int row, int column);
		unsigned char getLaplacian(int row, int column, ResponseLayer &src);
		float getResponse(int row, int column);
		float getResponse(int row, int column, ResponseLayer &src);
	};

	// These are passed in
	float thresh;
	int octaves;
	int init_sample;
	IntegralImage * img;

	// These get built
	std::vector<IPoint> * ipts;
	std::vector<ResponseLayer *> * responseMap;

	void buildResponseMap();
	void buildResponseLayer(ResponseLayer &rl);
	bool isExtremum(int r, int c, ResponseLayer &t, ResponseLayer &m, ResponseLayer &b);
	void interpolateExtremum(int r, int c, ResponseLayer &t, ResponseLayer &m, ResponseLayer &b);
	double * BuildDerivative(int r, int c, ResponseLayer &t, ResponseLayer &m, ResponseLayer &b);
	double * BuildHessian(int r, int c, ResponseLayer &t, ResponseLayer &m, ResponseLayer &b);
	double * Inverse(double * m);
	double * MMM_neg_3x3_3x1(double * A, double * B);
};

#endif
