#ifndef SURFDESCRIPTOR_H
#define SURFDESCRIPTOR_H

#include "IntegralImage.h"
#include "IPoint.h"
#include <vector>

class SurfDescriptor
{
private:
	// Gaussian distribution with sigma = 2.5.  Used as a fast lookup
	static float gauss25[49];
	float * d_gauss25;

	// The integral image which is being used
	IntegralImage * img;

public:
	SurfDescriptor() { }

	// Static one-call do it all function
	static void DecribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage *img);
	void DescribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage *img);
	void GetOrientation(IPoint &ip, float * d_resX, float * d_resY);
	void GetDescriptor(IPoint &ip, bool bUpright, bool bExtended, float * d_descriptor, float * d_len);
	double GetAngle(float X, float Y);
	float Gaussian(int x, int y, float sig);
	float Gaussian(float x, float y, float sig);
};

#endif
