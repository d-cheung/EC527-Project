#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include "IPoint.h"
#include "IntegralImage.h"
#include "FastHessian.h"


/// <summary>
/// Static one-call do it all function
/// </summary>
/// <param name="img"></param>
/// <param name="ipts"></param>
/// <param name="extended"></param>
/// <param name="upright"></param>
static void DecribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage &img)
{
    SurfDescriptor des = new SurfDescriptor();
    des->DescribeInterestPoints(ipts, upright, extended, img);
}


/// <summary>
    /// Build descriptor vector for each interest point in the supplied list
    /// </summary>
    /// <param name="img"></param>
    /// <param name="ipts"></param>
    /// <param name="upright"></param>
void DescribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage &img)
{
    if (ipts->size() == 0) return;
    this->img = img;
      
	for (std::vector<IPoint>::iterator ip = ipts->begin() ; ip != ipts->end(); ++ip)
    {
		// determine descriptor size
		if (extended) ip.descriptorLength = 128;
		else ip.descriptorLength = 64;

		// if we want rotation invariance get the orientation
		if (!upright) GetOrientation(ip);

		// Extract SURF descriptor
		GetDescriptor(ip, upright, extended);
    }
}

/// <summary>
    /// Determine dominant orientation for InterestPoint
    /// </summary>
    /// <param name="ip"></param>
void GetOrientation(IPoint &ip)
{
  const unsigned char Responses = 109;
  float[] resX = new float[Responses];
  float[] resY = new float[Responses];
  float[] Ang = new float[Responses];
  int idx = 0;
  int[] id = { 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6 };

  // Get rounded InterestPoint data
  int X = (int)floor(ip.x + (float)0.5);
  int Y = (int)floor(ip.y + (float)0.5);
  int S = (int)floor(ip.scale + (float)0.5);

  // calculate haar responses for points within radius of 6*scale
  for (int i = -6; i <= 6; ++i)
  {
	for (int j = -6; j <= 6; ++j)
	{
	  if (i * i + j * j < 36)
	  {
		float gauss = gauss25[id[i + 6], id[j + 6]];
		resX[idx] = gauss * img.HaarX(Y + j * S, X + i * S, 4 * S);
		resY[idx] = gauss * img.HaarY(Y + j * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	  }
	}
  }

  // calculate the dominant direction 
  float sumX, sumY, max = 0, orientation = 0;
  float ang1, ang2;
  float pi = (float)M_PI;

  // loop slides pi/3 window around feature point
  for (ang1 = 0; ang1 < 2 * pi; ang1 += (float)0.15)
  {
	ang2 = (ang1 + pi / (float)3 > 2 * pi ? ang1 - 5 * pi / (float)3 : ang1 + pi / (float)3);
	sumX = sumY = 0;

	for (int k = 0; k < Responses; ++k)
	{
	  // determine whether the point is within the window
	  if (ang1 < ang2 && ang1 < Ang[k] && Ang[k] < ang2)
	  {
		sumX += resX[k];
		sumY += resY[k];
	  }
	  else if (ang2 < ang1 &&
		((Ang[k] > 0 && Ang[k] < ang2) || (Ang[k] > ang1 && Ang[k] < pi)))
	  {
		sumX += resX[k];
		sumY += resY[k];
	  }
	}

	// if the vector produced from this window is longer than all 
	// previous vectors then this forms the new dominant direction
	if (sumX * sumX + sumY * sumY > max)
	{
	  // store largest orientation
	  max = sumX * sumX + sumY * sumY;
	  orientation = (float)GetAngle(sumX, sumY);
	}
  }

  // assign orientation of the dominant response vector
  ip.orientation = (float)orientation;
}


/// <summary>
    /// Construct descriptor vector for this interest point
    /// </summary>
    /// <param name="bUpright"></param>
void GetDescriptor(IPoint &ip, bool bUpright, bool bExtended)
{
    int sample_x, sample_y, count = 0;
    int i = 0, ix = 0, j = 0, jx = 0, xs = 0, ys = 0;
    float dx, dy, mdx, mdy, co, si;
    float dx_yn, mdx_yn, dy_xn, mdy_xn;
    float gauss_s1 = (float)0, gauss_s2 = (float)0;
    float rx = (float)0, ry = (float)0, rrx = (float)0, rry = (float)0, len = (float)0;
    float cx = (float)-0.5, cy = (float)0; //Subregion centers for the 4x4 gaussian weighting

    // Get rounded InterestPoint data
    int X = (int)floor(ip.x + (float)0.5);
    int Y = (int)floor(ip.y + (float)0.5);
    int S = (int)floor(ip.scale + (float)0.5);

    // Allocate descriptor memory
    ip.SetDescriptorLength(64);

    if (bUpright)
    {
		co = 1;
		si = 0;
    }
    else
    {
		co = (float)cos(ip.orientation);
		si = (float)sin(ip.orientation);
    }

    //Calculate descriptor for this interest point
    i = -8;
    while (i < 12)
    {
		j = -8;
		i = i - 4;

		cx += (float)1;
		cy = (float)-0.5;

		while (j < 12)
		{
			cy += (float)1;

			j = j - 4;

			ix = i + 5;
			jx = j + 5;

			dx = dy = mdx = mdy = (float)0;
			dx_yn = mdx_yn = dy_xn = mdy_xn = (float)0;

			xs = (int)floor(X + (-jx * S * si + ix * S * co) + (float)0.5);
			ys = (int)floor(Y + (jx * S * co + ix * S * si) + (float)0.5);

			// zero the responses
			dx = dy = mdx = mdy = (float)0;
			dx_yn = mdx_yn = dy_xn = mdy_xn = (float)0;

			for (int k = i; k < i + 9; ++k)
			{
				for (int l = j; l < j + 9; ++l)
				{
					//Get coords of sample point on the rotated axis
					sample_x = (int)floor(X + (-l * S * si + k * S * co) + (float)0.5);
					sample_y = (int)floor(Y + (l * S * co + k * S * si) + (float)0.5);

					//Get the gaussian weighted x and y responses
					gauss_s1 = Gaussian(xs - sample_x, ys - sample_y, (float)2.5 * S);
					rx = (float)img.HaarX(sample_y, sample_x, 2 * S);
					ry = (float)img.HaarY(sample_y, sample_x, 2 * S);

					//Get the gaussian weighted x and y responses on rotated axis
					rrx = gauss_s1 * (-rx * si + ry * co);
					rry = gauss_s1 * (rx * co + ry * si);


					if (bExtended)
					{
						// split x responses for different signs of y
						if (rry >= 0)
						{
							dx += rrx;
							mdx += abs(rrx);
						}
						else
						{
							dx_yn += rrx;
							mdx_yn += abs(rrx);
						}

						// split y responses for different signs of x
						if (rrx >= 0)
						{
							dy += rry;
							mdy += abs(rry);
						}
						else
						{
							dy_xn += rry;
							mdy_xn += abs(rry);
						}
					}
					else
					{
						dx += rrx;
						dy += rry;
						mdx += abs(rrx);
						mdy += abs(rry);
					}
				}
			}

			//Add the values to the descriptor vector
			gauss_s2 = Gaussian(cx - (float)2, cy - (float)2, (float)1.5);

			ip.descriptor[count++] = dx * gauss_s2;
			ip.descriptor[count++] = dy * gauss_s2;
			ip.descriptor[count++] = mdx * gauss_s2;
			ip.descriptor[count++] = mdy * gauss_s2;

			// add the extended components
			if (bExtended)
			{
				ip.descriptor[count++] = dx_yn * gauss_s2;
				ip.descriptor[count++] = dy_xn * gauss_s2;
				ip.descriptor[count++] = mdx_yn * gauss_s2;
				ip.descriptor[count++] = mdy_xn * gauss_s2;
			}

			len += (dx * dx + dy * dy + mdx * mdx + mdy * mdy
					+ dx_yn + dy_xn + mdx_yn + mdy_xn) * gauss_s2 * gauss_s2;

			j += 9;
		}
		i += 9;
    }

    //Convert to Unit Vector
    len = (float)sqrt((double)len);
    if (len > 0)
    {
		for (int d = 0; d < ip.descriptorLength; ++d)
		{
			ip.descriptor[d] /= len;
		}
    }
}


/// <summary>
/// Get the angle formed by the vector [x,y]
/// </summary>
/// <param name="X"></param>
/// <param name="Y"></param>
/// <returns></returns>
double GetAngle(float X, float Y)
{
  if (X >= 0 && Y >= 0)
	return atan(Y / X);
  else if (X < 0 && Y >= 0)
	return M_PI - atan(-Y / X);
  else if (X < 0 && Y < 0)
	return M_PI + atan(Y / X);
  else if (X >= 0 && Y < 0)
	return 2 * M_PI - atan(-Y / X);

  return 0;
}


/// <summary>
/// Get the value of the gaussian with std dev sigma
/// at the point (x,y)
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="sig"></param>
/// <returns></returns>
float Gaussian(int x, int y, float sig)
{
  float pi = (float)M_PI;
  return ((float)1 / ((float)2 * pi * sig * sig)) * (float)exp(-(x * x + y * y) / ((float)2.0 * sig * sig));
}


/// <summary>
/// Get the value of the gaussian with std dev sigma
/// at the point (x,y)
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="sig"></param>
/// <returns></returns>
float Gaussian(float x, float y, float sig)
{
  float pi = (float)M_PI;
  return (float)1 / ((float)2 * pi * sig * sig) * (float)exp(-(x * x + y * y) / ((float)2.0 * sig * sig));
}

