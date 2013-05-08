#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include "SurfDescriptor.h"
#include "IPoint.h"
#include "IntegralImage.h"
#include "FastHessian.h"
#include "parameters.h"
#include <pthread.h>

struct thread_data {
	int thread_id;
	IntegralImage* img;
	int i,j, xs, ys;
	int X,Y,S;
	float dx, dy, mdx, mdy, co, si;
    	float dx_yn, mdx_yn, dy_xn, mdy_xn;
	bool bExtended;
};

float SurfDescriptor::gauss25[7][7] = {
	{(float)0.02350693969273,(float)0.01849121369071,(float)0.01239503121241,(float)0.00708015417522,(float)0.00344628101733,(float)0.00142945847484,(float)0.00050524879060}, 
	{(float)0.02169964028389,(float)0.01706954162243,(float)0.01144205592615,(float)0.00653580605408,(float)0.00318131834134,(float)0.00131955648461,(float)0.00046640341759}, 
	{(float)0.01706954162243,(float)0.01342737701584,(float)0.00900063997939,(float)0.00514124713667,(float)0.00250251364222,(float)0.00103799989504,(float)0.00036688592278},
	{(float)0.01144205592615,(float)0.00900063997939,(float)0.00603330940534,(float)0.00344628101733,(float)0.00167748505986,(float)0.00069579213743,(float)0.00024593098864},
	{(float)0.00653580605408,(float)0.00514124713667,(float)0.00344628101733,(float)0.00196854695367,(float)0.00095819467066,(float)0.00039744277546,(float)0.00014047800980},
	{(float)0.00318131834134,(float)0.00250251364222,(float)0.00167748505986,(float)0.00095819467066,(float)0.00046640341759,(float)0.00019345616757,(float)0.00006837798818},
	{(float)0.00131955648461,(float)0.00103799989504,(float)0.00069579213743,(float)0.00039744277546,(float)0.00019345616757,(float)0.00008024231247,(float)0.00002836202103}
};


// Static one-call do it all function
void SurfDescriptor::DecribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage * img)
{
    SurfDescriptor des;
    des.DescribeInterestPoints(ipts, upright, extended, img);
}


// Build descriptor vector for each interest point in the supplied list
void SurfDescriptor::DescribeInterestPoints(std::vector<IPoint>* ipts, bool upright, bool extended, IntegralImage *img)
{
    if (ipts->size() == 0) return;
    this->img = img;
      
	for (std::vector<IPoint>::iterator ip = ipts->begin(); ip != ipts->end(); ++ip)
    {
		// determine descriptor size
		if (extended) ip->descriptorLength = 128;
		else ip->descriptorLength = 64;

		// if we want rotation invariance get the orientation
		if (!upright) GetOrientation(*ip);

		// Extract SURF descriptor
		GetDescriptor(*ip, upright, extended);
    }
}

// Determine dominant orientation for InterestPoint
void SurfDescriptor::GetOrientation(IPoint &ip)
{
  const unsigned char Responses = 109;
  float resX[Responses];
  float resY[Responses];
  float Ang[Responses];
  int idx = 0;
  int id[13] = { 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6 };

  // Get rounded InterestPoint data
  int X = (int)floor(ip.x + (float)0.5);
  int Y = (int)floor(ip.y + (float)0.5);
  int S = (int)floor(ip.scale + (float)0.5);

  // calculate haar responses for points within radius of 6*scale
  for (int i = -6; i <= 6; ++i)
  {
	int j= -6;
	int ii = i*i;
	int j1 = (j+0)*(j+0);
	int j2 = (j+1)*(j+1);
	int j3 = (j+2)*(j+2);
	int j4 = (j+3)*(j+3);
	int j5 = (j+4)*(j+4);
	int j6 = (j+5)*(j+5);
	int j7 = (j+6)*(j+6);
	int j8 = (j+7)*(j+7);
	int j9 = (j+8)*(j+8);
	int j10= (j+9)*(j+9);
	int j11= (j+10)*(j+10);
	int j12= (j+11)*(j+11);
	int j13= (j+12)*(j+12);


	if (ii + j1 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+0 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+0) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+0) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j2 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+1 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+1) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+1) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j3 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+2 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+2) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+2) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j4 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+3 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+3) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+3) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j5 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+4 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+4) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+4) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j6 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+5 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+5) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+5) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j7 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+6 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+6) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+6) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j8 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+7 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+7) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+7) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j9 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+8 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+8) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+8) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j10 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+9 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+9) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+9) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j11 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+10 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+10) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+10) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j12 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+11 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+11) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+11) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
	}
	if (ii + j13 < 36)
	{
		float gauss = gauss25[id[i + 6]][id[j+12 + 6]];
		resX[idx] = gauss * img->HaarX(Y + (j+12) * S, X + i * S, 4 * S);
		resY[idx] = gauss * img->HaarY(Y + (j+12) * S, X + i * S, 4 * S);
		Ang[idx] = (float)GetAngle(resX[idx], resY[idx]);
		++idx;
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


// Construct descriptor vector for this interest point
void SurfDescriptor::GetDescriptor(IPoint &ip, bool bUpright, bool bExtended)
{
    void * GetDescriptor_work(void * threadarg);
    pthread_t threads[NUM_THREADS];
    struct thread_data data[NUM_THREADS];

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
    for(int i = -12; i < 4; i+=5)
    //while (i < 12)
    {
		j = -8;
		//i = i - 4;

		cx += (float)1;
		cy = (float)-0.5;
		
		for(int j = -12; j < 4; j+=5)
		//while (j < 12)
		{
			cy += (float)1;

			//j = j - 4;
			//printf("i=%d, j=%d\n",i,j);
			ix = i + 5;
			jx = j + 5;

			xs = (int)floor(X + (-jx * S * si + ix * S * co) + (float)0.5);
			ys = (int)floor(Y + (jx * S * co + ix * S * si) + (float)0.5);

			// zero the responses
			dx = dy = mdx = mdy = (float)0;
			dx_yn = mdx_yn = dy_xn = mdy_xn = (float)0;


			for (int t = 0; t <NUM_THREADS; t++)
			{
				data[t].thread_id = t;
				data[t].img = img;
				data[t].i = i;
				data[t].j = j;
				data[t].S = S;
				data[t].co = co;
				data[t].si = si;
				data[t].X = X;
				data[t].Y = Y;
				data[t].dx = 0;
				data[t].mdx =0;
				data[t].dx_yn = 0;
				data[t].mdx_yn = 0;
				data[t].dy = 0;
				data[t].mdy = 0;
				data[t].dy_xn = 0;
				data[t].mdy_xn = 0;
				data[t].xs = xs;
				data[t].ys = ys;
				data[t].bExtended = bExtended;

				int rc = pthread_create(&threads[t], NULL, GetDescriptor_work, (void*) &(data[t]));
				if (rc)
				{
					printf("ERROR: return code from pthread_create() is %d\n", rc);
					exit(-1);
				}
			}
			for (int t = 0; t < NUM_THREADS; t++)
			{
				if (pthread_join(threads[t], NULL))
				{
					printf("\n ERROR on join SURFDESCRIPTOR\n");
					exit(19);
				}
				dx += data[t].dx;
				mdx += data[t].mdx;
				dx_yn += data[t].dx_yn;
				mdx_yn += data[t].mdx_yn;
				dy += data[t].dy;
				mdy += data[t].mdy;
				dy_xn += data[t].dy_xn;
				mdy_xn += data[t].mdy_xn;				
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

			//j += 9;
		}
		//i += 9;
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


void * GetDescriptor_work(void* threadarg)
{
	float myGaussian(int x, int y, float sig);
	float myHaarX(int row, int column, int size, IntegralImage* img);
	float myHaarY(int row, int column, int size, IntegralImage* img);
	thread_data* data = (thread_data *) threadarg;
	int t = data->thread_id;
	IntegralImage* img = data->img;
	int i = data->i;
	int j = data->j;
	int S = data->S;
	float co = data->co;
	float si = data->si;
	int X = data->X;
	int Y = data->Y;
	int xs = data->xs;
	int ys = data->ys;
	bool bExtended = data->bExtended;
	
	int start,end;
	if (NUM_THREADS <9)
	{
		start = NUM_THREADS * t + i;
		end = NUM_THREADS + start;
	}
	else
	{
		start = t + i;
		end = 1 + start;
	}
	if (end > i + 9) end  = i + 9;

	//printf("Thread %d, i %d, start %d, end %d\n",t,i, start,end);
	for (int k = start; k < end; ++k)
	{
		int l;
				for (l = j; l < j + 9 - 4; l+=4)
				{
					//Get coords of sample point on the rotated axis
					int sample_x0 = (int)floor(X + ((-l+0) * S * si + k * S * co) + (float)0.5);
					int sample_x1 = (int)floor(X + ((-l+1) * S * si + k * S * co) + (float)0.5);
					int sample_x2 = (int)floor(X + ((-l+2) * S * si + k * S * co) + (float)0.5);
					int sample_x3 = (int)floor(X + ((-l+3) * S * si + k * S * co) + (float)0.5);

					int sample_y0 = (int)floor(Y + ((l+0) * S * co + k * S * si) + (float)0.5);
					int sample_y1 = (int)floor(Y + ((l+1) * S * co + k * S * si) + (float)0.5);
					int sample_y2 = (int)floor(Y + ((l+2) * S * co + k * S * si) + (float)0.5);
					int sample_y3 = (int)floor(Y + ((l+3) * S * co + k * S * si) + (float)0.5);

					//Get the gaussian weighted x and y responses
					float gauss_s1_0 = myGaussian(xs - sample_x0, ys - sample_y0, (float)2.5 * S);
					float gauss_s1_1 = myGaussian(xs - sample_x1, ys - sample_y1, (float)2.5 * S);
					float gauss_s1_2 = myGaussian(xs - sample_x2, ys - sample_y2, (float)2.5 * S);
					float gauss_s1_3 = myGaussian(xs - sample_x3, ys - sample_y3, (float)2.5 * S);

					float rx0 = (float)img->HaarX(sample_y0, sample_x0, 2 * S);
					float rx1 = (float)img->HaarX(sample_y1, sample_x1, 2 * S);
					float rx2 = (float)img->HaarX(sample_y2, sample_x2, 2 * S);
					float rx3 = (float)img->HaarX(sample_y3, sample_x3, 2 * S);

					float ry0 = (float)img->HaarY(sample_y0, sample_x0, 2 * S);
					float ry1 = (float)img->HaarY(sample_y1, sample_x1, 2 * S);
					float ry2 = (float)img->HaarY(sample_y2, sample_x2, 2 * S);
					float ry3 = (float)img->HaarY(sample_y3, sample_x3, 2 * S);

					//Get the gaussian weighted x and y responses on rotated axis
					float rrx0 = gauss_s1_0 * (-rx0 * si + ry0 * co);
					float rrx1 = gauss_s1_1 * (-rx1 * si + ry1 * co);
					float rrx2 = gauss_s1_2 * (-rx2 * si + ry2 * co);
					float rrx3 = gauss_s1_3 * (-rx3 * si + ry3 * co);

					float rry0 = gauss_s1_0 * (rx0 * co + ry0 * si);
					float rry1 = gauss_s1_1 * (rx1 * co + ry1 * si);
					float rry2 = gauss_s1_2 * (rx2 * co + ry2 * si);
					float rry3 = gauss_s1_3 * (rx3 * co + ry3 * si);

					if (bExtended)
					{
						// split x responses for different signs of y
						if (rry0 >= 0)
						{
							data->dx += rrx0;
							data->mdx += fabs(rrx0);
						}
						else
						{
							data->dx_yn += rrx0;
							data->mdx_yn += fabs(rrx0);
						}

						// split y responses for different signs of x
						if (rrx0 >= 0)
						{
							data->dy += rry0;
							data->mdy += fabs(rry0);
						}
						else
						{
							data->dy_xn += rry0;
							data->mdy_xn += fabs(rry0);
						}
						// split x responses for different signs of y
						if (rry1 >= 0)
						{
							data->dx += rrx1;
							data->mdx += fabs(rrx1);
						}
						else
						{
							data->dx_yn += rrx1;
							data->mdx_yn += fabs(rrx1);
						}

						// split y responses for different signs of x
						if (rrx1 >= 0)
						{
							data->dy += rry1;
							data->mdy += fabs(rry1);
						}
						else
						{
							data->dy_xn += rry1;
							data->mdy_xn += fabs(rry1);
						}
						// split x responses for different signs of y
						if (rry2 >= 0)
						{
							data->dx += rrx2;
							data->mdx += fabs(rrx2);
						}
						else
						{
							data->dx_yn += rrx2;
							data->mdx_yn += fabs(rrx2);
						}

						// split y responses for different signs of x
						if (rrx2 >= 0)
						{
							data->dy += rry2;
							data->mdy += fabs(rry2);
						}
						else
						{
							data->dy_xn += rry2;
							data->mdy_xn += fabs(rry2);
						}
						// split x responses for different signs of y
						if (rry3 >= 0)
						{
							data->dx += rrx3;
							data->mdx += fabs(rrx3);
						}
						else
						{
							data->dx_yn += rrx3;
							data->mdx_yn += fabs(rrx3);
						}

						// split y responses for different signs of x
						if (rrx3 >= 0)
						{
							data->dy += rry3;
							data->mdy += fabs(rry3);
						}
						else
						{
							data->dy_xn += rry3;
							data->mdy_xn += fabs(rry3);
						}
					}
					else
					{
						data->dx += rrx0;
						data->dx += rrx1;
						data->dx += rrx2;
						data->dx += rrx3;

						data->dy += rry0;
						data->dy += rry1;
						data->dy += rry2;
						data->dy += rry3;

						data->mdx += fabs(rrx0);
						data->mdx += fabs(rrx1);
						data->mdx += fabs(rrx2);
						data->mdx += fabs(rrx3);

						data->mdy += fabs(rry0);
						data->mdy += fabs(rry1);
						data->mdy += fabs(rry2);
						data->mdy += fabs(rry3);
					}
				}
				for (; l < j + 9; ++l)
				{
					//Get coords of sample point on the rotated axis
					int sample_x0 = (int)floor(X + ((-l+0) * S * si + k * S * co) + (float)0.5);
					int sample_y0 = (int)floor(Y + ((l+0) * S * co + k * S * si) + (float)0.5);

					//Get the gaussian weighted x and y responses
					float gauss_s1_0 = myGaussian(xs - sample_x0, ys - sample_y0, (float)2.5 * S);
					float rx0 = (float)img->HaarX(sample_y0, sample_x0, 2 * S);
					float ry0 = (float)img->HaarY(sample_y0, sample_x0, 2 * S);

					//Get the gaussian weighted x and y responses on rotated axis
					float rrx0 = gauss_s1_0 * (-rx0 * si + ry0 * co);
					float rry0 = gauss_s1_0 * (rx0 * co + ry0 * si);

					if (bExtended)
					{
						// split x responses for different signs of y
						if (rry0 >= 0)
						{
							data->dx += rrx0;
							data->mdx += fabs(rrx0);
						}
						else
						{
							data->dx_yn += rrx0;
							data->mdx_yn += fabs(rrx0);
						}

						// split y responses for different signs of x
						if (rrx0 >= 0)
						{
							data->dy += rry0;
							data->mdy += fabs(rry0);
						}
						else
						{
							data->dy_xn += rry0;
							data->mdy_xn += fabs(rry0);
						}
					}
					else
					{
						data->dx += rrx0;
						data->dy += rry0;
						data->mdx += fabs(rrx0);
						data->mdy += fabs(rry0);
					}
				}
	}
	pthread_exit(NULL);
}

float myGaussian(int x, int y, float sig)
{
  float pi = (float)M_PI;
  return ((float)1 / ((float)2 * pi * sig * sig)) * (float)exp(-(x * x + y * y) / ((float)2.0 * sig * sig));
}


// Get the angle formed by the vector [x,y]
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


// Get the value of the gaussian with std dev sigma
// at the point (x,y)
float SurfDescriptor::Gaussian(int x, int y, float sig)
{
  float pi = (float)M_PI;
  return ((float)1 / ((float)2 * pi * sig * sig)) * (float)exp(-(x * x + y * y) / ((float)2.0 * sig * sig));
}


// Get the value of the gaussian with std dev sigma
// at the point (x,y)
float SurfDescriptor::Gaussian(float x, float y, float sig)
{
  float pi = (float)M_PI;
  return (float)1 / ((float)2 * pi * sig * sig) * (float)exp(-(x * x + y * y) / ((float)2.0 * sig * sig));
}

