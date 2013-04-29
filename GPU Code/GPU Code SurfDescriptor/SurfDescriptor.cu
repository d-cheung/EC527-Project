#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include "SurfDescriptor.h"
#include "IPoint.h"
#include "IntegralImage.h"
#include "FastHessian.h"

float SurfDescriptor::gauss25[49] = {
	(float)0.02350693969273,(float)0.01849121369071,(float)0.01239503121241,(float)0.00708015417522,(float)0.00344628101733,(float)0.00142945847484,(float)0.00050524879060, 
	(float)0.02169964028389,(float)0.01706954162243,(float)0.01144205592615,(float)0.00653580605408,(float)0.00318131834134,(float)0.00131955648461,(float)0.00046640341759, 
	(float)0.01706954162243,(float)0.01342737701584,(float)0.00900063997939,(float)0.00514124713667,(float)0.00250251364222,(float)0.00103799989504,(float)0.00036688592278,
	(float)0.01144205592615,(float)0.00900063997939,(float)0.00603330940534,(float)0.00344628101733,(float)0.00167748505986,(float)0.00069579213743,(float)0.00024593098864,
	(float)0.00653580605408,(float)0.00514124713667,(float)0.00344628101733,(float)0.00196854695367,(float)0.00095819467066,(float)0.00039744277546,(float)0.00014047800980,
	(float)0.00318131834134,(float)0.00250251364222,(float)0.00167748505986,(float)0.00095819467066,(float)0.00046640341759,(float)0.00019345616757,(float)0.00006837798818,
	(float)0.00131955648461,(float)0.00103799989504,(float)0.00069579213743,(float)0.00039744277546,(float)0.00019345616757,(float)0.00008024231247,(float)0.00002836202103
};


/// <summary>
/// Static one-call do it all function
/// </summary>
/// <param name="img"></param>
/// <param name="ipts"></param>
/// <param name="extended"></param>
/// <param name="upright"></param>
void SurfDescriptor::DecribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage * img)
{
    SurfDescriptor des;
    des.DescribeInterestPoints(ipts, upright, extended, img);
}


/// <summary>
    /// Build descriptor vector for each interest point in the supplied list
    /// </summary>
    /// <param name="img"></param>
    /// <param name="ipts"></param>
    /// <param name="upright"></param>
void SurfDescriptor::DescribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage *img)
{
	if (ipts->size() == 0) return;
	this->img = img;

	float * d_resX, * d_resY, * d_Ang;
	float * d_descriptor, * d_len;

	cudaMalloc((void**) &d_resX,      (109)*sizeof(float));
	cudaMalloc((void**) &d_resY,      (109)*sizeof(float));
	cudaMalloc((void**) &d_Ang,       (109)*sizeof(float));
	cudaMalloc((void**) &d_gauss25,    (49)*sizeof(float));
	cudaMalloc((void**) &d_descriptor, (64)*sizeof(float));
	cudaMalloc((void**) &d_len,        (16)*sizeof(float));

	cudaMemcpy(d_gauss25, gauss25, 49*sizeof(float), cudaMemcpyHostToDevice);

	for (std::vector<IPoint>::iterator ip = ipts->begin(); ip != ipts->end(); ++ip)
	{
		// determine descriptor size
		if (extended) ip->descriptorLength = 128;
		else ip->descriptorLength = 64;

		// if we want rotation invariance get the orientation
		if (!upright) GetOrientation(*ip, d_resX, d_resY, d_Ang);

		// Extract SURF descriptor
		GetDescriptor(*ip, upright, extended, d_descriptor, d_len);
	}

	cudaFree(d_resX);
	cudaFree(d_resY);
	cudaFree(d_Ang);
	cudaFree(d_gauss25);
	cudaFree(d_len);
	cudaFree(d_descriptor);
}
/*
__device__ double d_GetAngle(float X, float Y)
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
*/

__global__ void cudaHaar6x6(float ** img, float * resX, float * resY, float * Ang, float * gauss25, int X, int Y, int S, int Height, int Width)
{
	int id[11] = {5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5};
	int height = Height;
	int width = Width;
	
	int i = id[threadIdx.x];
	int j = id[threadIdx.y];

	int idx;
	switch (threadIdx.x)
	{
		case 0:
			switch (threadIdx.y) {
				case 0:
				case 1:
				case 9:
				case 10:
					return;
				default:
					idx = threadIdx.y - 2;
			}
			break;
		case 1:
			switch (threadIdx.y) {
				case 0:
				case 10:
					return;
				default:
					idx = threadIdx.y + 6;
			}
			break;
		case 9:
			switch (threadIdx.y) {
				case 0:
				case 10:
					return;
				default:
					idx = threadIdx.y + 92;
			}
			break;
		case 10:
			switch (threadIdx.y) {
				case 0:
				case 1:
				case 9:
				case 10:
					return;
				default:
					idx = threadIdx.y + 101;
			}
			break;
		default:
			idx = threadIdx.x * 11 + threadIdx.y - 6;
	}

	float gauss = gauss25[id[i]*7 + id[j]];
	resX[idx] = gauss * d_HaarX(img, Y + j * S, X + i * S, 4 * S, height, width);
	resY[idx] = gauss * d_HaarY(img, Y + j * S, X + i * S, 4 * S, height, width);
//	Ang[idx] = (float)d_GetAngle(resX[idx], resY[idx]);
}

__global__ void cudaHaar6x6(float * img, float * resX, float * resY, float * Ang, float * gauss25, int X, int Y, int S, size_t Pitch, int Height, int Width)
{
	int id[11] = {5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5};
	int pitch = Pitch;
	int height = Height;
	int width = Width;
	
	int i = id[threadIdx.x];
	int j = id[threadIdx.y];

	int idx;
	switch (threadIdx.x)
	{
		case 0:
			switch (threadIdx.y) {
				case 0:
				case 1:
				case 9:
				case 10:
					return;
				default:
					idx = threadIdx.y - 2;
			}
			break;
		case 1:
			switch (threadIdx.y) {
				case 0:
				case 10:
					return;
				default:
					idx = threadIdx.y + 6;
			}
			break;
		case 9:
			switch (threadIdx.y) {
				case 0:
				case 10:
					return;
				default:
					idx = threadIdx.y + 92;
			}
			break;
		case 10:
			switch (threadIdx.y) {
				case 0:
				case 1:
				case 9:
				case 10:
					return;
				default:
					idx = threadIdx.y + 101;
			}
			break;
		default:
			idx = threadIdx.x * 11 + threadIdx.y - 6;
	}

	float gauss = gauss25[id[i]*7 + id[j]];
	resX[idx] = gauss * d_HaarX(img, Y + j * S, X + i * S, 4 * S, pitch, height, width);
	resY[idx] = gauss * d_HaarY(img, Y + j * S, X + i * S, 4 * S, pitch, height, width);
//	Ang[idx] = (float)d_GetAngle(resX[idx], resY[idx]);
  }

/// <summary>
    /// Determine dominant orientation for InterestPoint
    /// </summary>
    /// <param name="ip"></param>
void SurfDescriptor::GetOrientation(IPoint &ip, float * d_resX, float * d_resY, float * d_Ang)
{
  const unsigned char Responses = 109;
  float resX[Responses];
  float resY[Responses];
  float Ang[Responses];

  // Get rounded InterestPoint data
  int X = (int)floor(ip.x + (float)0.5);
  int Y = (int)floor(ip.y + (float)0.5);
  int S = (int)floor(ip.scale + (float)0.5);

  dim3 dimBlock(11, 11);


#ifdef FLOAT_SINGLE
  cudaHaar6x6 <<<1, dimBlock>>> (img->Matrix, d_resX, d_resY, d_Ang, d_gauss25, X, Y, S, img->Pitch, img->Height, img->Width);
#endif

#ifdef FLOAT_DOUBLE
  cudaHaar6x6 <<<1, dimBlock>>> (img->Matrix, d_resX, d_resY, d_Ang, d_gauss25, X, Y, S, img->Height, img->Width);
#endif

  cudaThreadSynchronize();

  cudaMemcpy(resX, d_resX, (Responses)*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(resY, d_resY, (Responses)*sizeof(float), cudaMemcpyDeviceToHost);
//  cudaMemcpy( Ang, d_Ang,  (Responses)*sizeof(float), cudaMemcpyDeviceToHost);

  for (unsigned char ii = 0; ii < Responses; ii++)
	Ang[ii] = (float)GetAngle(resX[ii], resY[ii]);

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

__device__ float d_Gaussian(int x, int y, float sig)
{
  float pi = (float)M_PI;
  return ((float)1 / ((float)2 * pi * sig * sig)) * (float)exp(-(x * x + y * y) / ((float)2.0 * sig * sig));
}

__global__ void cudaGetDescriptor(float ** img, float * descriptor, float * len, int d_X, int d_Y, int d_S, float d_co, float d_si, int Height, int Width)
{
	extern __shared__ float floatdata[];

	int X = d_X;
	int Y = d_Y;
	int S = d_S;
	float co = d_co;
	float si = d_si;

	int ij_id[4] = {-12, -7, -2, 3};

	int i = ij_id[blockIdx.x];
	int j = ij_id[blockIdx.y];
	int bid = blockIdx.x*gridDim.x + blockIdx.y;
	float cx = (float)(-0.5) + (float)(1 + blockIdx.x);
	float cy = (float)(-0.5) + (float)(1 + blockIdx.y);

	int ix = i + 5;
	int jx = j + 5;

	float dx, dy, mdx, mdy;

	float xs = (int)floor(X + (-jx * S * si + ix * S * co) + (float)0.5);
	float ys = (int)floor(Y + (jx * S * co + ix * S * si) + (float)0.5);


	int k = i + threadIdx.x;
	int l = j + threadIdx.y;
	int height = Height;
	int width = Width;
	int idx = threadIdx.x*blockDim.x + threadIdx.y;

	//Get coords of sample point on the rotated axis
	int sample_x = (int)floor(X + (-l * S * si + k * S * co) + (float)0.5);
	int sample_y = (int)floor(Y + (l * S * co + k * S * si) + (float)0.5);

	//Get the gaussian weighted x and y responses
	float gauss_s1 = d_Gaussian(xs - sample_x, ys - sample_y, (float)2.5 * S);
	float rx = (float)d_HaarX(img, sample_y, sample_x, 2 * S, height, width);
	float ry = (float)d_HaarY(img, sample_y, sample_x, 2 * S, height, width);

	//Get the gaussian weighted x and y responses on rotated axis
	float rrx = gauss_s1 * (-rx * si + ry * co);
	float rry = gauss_s1 * (rx * co + ry * si);

	floatdata[0*81+idx] = rrx;
	floatdata[1*81+idx] = rry;
	floatdata[2*81+idx] = fabs(rrx);
	floatdata[3*81+idx] = fabs(rry);

	__syncthreads();

	if (idx == 39)
	{
		floatdata[0*81+39] += floatdata[0*81+80] + floatdata[0*81+79];
		floatdata[1*81+39] += floatdata[1*81+80] + floatdata[1*81+79];
		floatdata[2*81+39] += floatdata[2*81+80] + floatdata[2*81+79];
		floatdata[3*81+39] += floatdata[3*81+80] + floatdata[3*81+79];
	}

	else if (idx < 39)
	{
		floatdata[0*81+idx] += floatdata[0*81+idx + 40];
		floatdata[1*81+idx] += floatdata[1*81+idx + 40];
		floatdata[2*81+idx] += floatdata[2*81+idx + 40];
		floatdata[3*81+idx] += floatdata[3*81+idx + 40];
	}

	__syncthreads();

	if (idx < 20)
	{
		floatdata[0*81+idx] += floatdata[0*81+idx + 20];
		floatdata[1*81+idx] += floatdata[1*81+idx + 20];
		floatdata[2*81+idx] += floatdata[2*81+idx + 20];
		floatdata[3*81+idx] += floatdata[3*81+idx + 20];
	}

	__syncthreads();

	if (idx < 10)
	{
		floatdata[0*81+idx] += floatdata[0*81+idx + 10];
		floatdata[1*81+idx] += floatdata[1*81+idx + 10];
		floatdata[2*81+idx] += floatdata[2*81+idx + 10];
		floatdata[3*81+idx] += floatdata[3*81+idx + 10];
	}

	__syncthreads();

	if (idx < 5)
	{
		floatdata[0*81+idx] += floatdata[0*81+idx + 5];
		floatdata[1*81+idx] += floatdata[1*81+idx + 5];
		floatdata[2*81+idx] += floatdata[2*81+idx + 5];
		floatdata[3*81+idx] += floatdata[3*81+idx + 5];
	}

	__syncthreads();

	if (idx == 0)
	{
		floatdata[0*81+0] += floatdata[0*81+1] + floatdata[0*81+2] + floatdata[0*81+3] + floatdata[0*81+4];
		floatdata[1*81+0] += floatdata[1*81+1] + floatdata[1*81+2] + floatdata[1*81+3] + floatdata[1*81+4];
		floatdata[2*81+0] += floatdata[2*81+1] + floatdata[2*81+2] + floatdata[2*81+3] + floatdata[2*81+4];
		floatdata[3*81+0] += floatdata[3*81+1] + floatdata[3*81+2] + floatdata[3*81+3] + floatdata[3*81+4];

		dx =  floatdata[0*81+0];
		dy =  floatdata[1*81+0];
		mdx = floatdata[2*81+0];
		mdy = floatdata[3*81+0];
		//Add the values to the descriptor vector
		float gauss_s2 = d_Gaussian(cx - (float)2, cy - (float)2, (float)1.5);

		descriptor[0+bid] =  dx * gauss_s2;
		descriptor[1+bid] =  dy * gauss_s2;
		descriptor[2+bid] = mdx * gauss_s2;
		descriptor[3+bid] = mdy * gauss_s2;
		len[bid] = (dx * dx + dy * dy + mdx * mdx + mdy * mdy) * gauss_s2 * gauss_s2;
	}
}


/// <summary>
    /// Construct descriptor vector for this interest point
    /// </summary>
    /// <param name="bUpright"></param>
void SurfDescriptor::GetDescriptor(IPoint &ip, bool bUpright, bool bExtended, float * d_descriptor, float * d_len)
{
    float co, si;
    float len = (float)0;

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


	dim3 dimGrid(4, 4);
	dim3 dimBlock(9, 9);


	cudaGetDescriptor <<<dimGrid, dimBlock, 4*81*sizeof(float)>>> (img->Matrix, d_descriptor, d_len, X, Y, S, co, si, img->Height, img->Width);
	cudaThreadSynchronize();


	float h_len[16];
	cudaMemcpy(ip.descriptor, d_descriptor, 64*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(        h_len,        d_len, 16*sizeof(float), cudaMemcpyDeviceToHost);


	for (int ii = 0; ii < 16; ii++)
		len += h_len[ii];

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
double SurfDescriptor::GetAngle(float X, float Y)
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
float SurfDescriptor::Gaussian(int x, int y, float sig)
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
float SurfDescriptor::Gaussian(float x, float y, float sig)
{
  float pi = (float)M_PI;
  return (float)1 / ((float)2 * pi * sig * sig) * (float)exp(-(x * x + y * y) / ((float)2.0 * sig * sig));
}

