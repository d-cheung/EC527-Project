#ifndef SURFDESCRIPTOR_H
#define SURFDESCRIPTOR_H

class SurfDescriptor
{
private:
	/// Gaussian distribution with sigma = 2.5.  Used as a fast lookup
	float[][] gauss25 = new float[7][7] {
		{0.02350693969273f,0.01849121369071f,0.01239503121241f,0.00708015417522f,0.00344628101733f,0.00142945847484f,0.00050524879060f},
		{0.02169964028389f,0.01706954162243f,0.01144205592615f,0.00653580605408f,0.00318131834134f,0.00131955648461f,0.00046640341759f},
		{0.01706954162243f,0.01342737701584f,0.00900063997939f,0.00514124713667f,0.00250251364222f,0.00103799989504f,0.00036688592278f},
		{0.01144205592615f,0.00900063997939f,0.00603330940534f,0.00344628101733f,0.00167748505986f,0.00069579213743f,0.00024593098864f},
		{0.00653580605408f,0.00514124713667f,0.00344628101733f,0.00196854695367f,0.00095819467066f,0.00039744277546f,0.00014047800980f},
		{0.00318131834134f,0.00250251364222f,0.00167748505986f,0.00095819467066f,0.00046640341759f,0.00019345616757f,0.00006837798818f},
		{0.00131955648461f,0.00103799989504f,0.00069579213743f,0.00039744277546f,0.00019345616757f,0.00008024231247f,0.00002836202103f}
	};
	/// The integral image which is being used
	IntegralImage img;

public:
	SurfDescriptor()
	{

	}
	static void DecribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage &img);
	void DescribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage &img);
	void GetOrientation(IPoint &ip);
	void GetDescriptor(IPoint &ip, bool bUpright, bool bExtended);
	double GetAngle(float X, float Y);
	float Gaussian(int x, int y, float sig);
	float Gaussian(float x, float y, float sig);
}; // SurfDescriptor


#endif