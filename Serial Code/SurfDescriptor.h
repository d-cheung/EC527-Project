#ifndef SURFDESCRIPTOR_H
#define SURFDESCRIPTOR_H

#include "IntegralImage.h"
#include "IPoint.h"
#include <vector>

class SurfDescriptor
{

private:
	// Gaussian distribution with sigma = 2.5.  Used as a fast lookup
	static float gauss25[7][7];
	// The integral image which is being used
	IntegralImage * img;

public:
	SurfDescriptor() { }
	static void DecribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage *img);
	void DescribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage *img);
	void GetOrientation(IPoint &ip);
	void GetDescriptor(IPoint &ip, bool bUpright, bool bExtended);
	double GetAngle(float X, float Y);
	float Gaussian(int x, int y, float sig);
	float Gaussian(float x, float y, float sig);
};

#endif
