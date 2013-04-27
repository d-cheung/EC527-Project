#include <iostream>
#include "FastHessian.h"
#include "IntegralImage.h"
#include "bitmap_image.hpp"
#include "IPoint.h"
#include "SurfDescriptor.h"
#include <vector>
#include <cmath>
#include "timespec.h"

#define CPG 2.53           // Cycles per GHz -- Adjust to your computer
#define ITERATIONS 1

/************************************************************************************
*  Usage:
*      ./surf [source file name] [destination file name]
*
************************************************************************************/

int main(int argc, char ** argv)
{
	if (argc == 3)
	{

		void PaintSURF(bitmap_image &img, std::vector<IPoint> * ipts);
		int clock_gettime(clockid_t clk_id, struct timespec *tp);

		struct timespec time1, time2;
		struct timespec time_stamp[3][ITERATIONS];
		long int i, j;
		const char * labels[3] = {"IntegralImage","FastHessian","SurfDescriptor"};

		std::string file_name(argv[1]);
		bitmap_image image(file_name);

		IntegralImage * iimg;
		// Running IntegralImage
		for (i = 0; i < ITERATIONS; i++)
		{
			clock_gettime(CLOCK_REALTIME, &time1);
			iimg = IntegralImage::FromImage(image);
			clock_gettime(CLOCK_REALTIME, &time2);
			time_stamp[0][i] = diff(time1,time2);
		}

		std::vector<IPoint> * ipts = NULL;

		// Running FastHessian
		for (i = 0; i < ITERATIONS; i++)
		{
			clock_gettime(CLOCK_REALTIME, &time1);
			ipts = FastHessian::getIpoints((float)0.0002, 5, 2, iimg);
			clock_gettime(CLOCK_REALTIME, &time2);
			time_stamp[1][i] = diff(time1,time2);
		}

		if (ipts == NULL)
		{
			printf("getIpoints failed\n");
			return 0;
		}

		// Running SurfDescriptor
		for (i = 0; i < ITERATIONS; i++)
		{
			clock_gettime(CLOCK_REALTIME, &time1);
			SurfDescriptor::DecribeInterestPoints(ipts, false, false, iimg);
			clock_gettime(CLOCK_REALTIME, &time2);
			time_stamp[2][i] = diff(time1,time2);
		}

		PaintSURF(image, ipts);
		image.save_image(argv[2]);

		printf("Width: %d, Height, %d\n", image.width(), image.height());
		for (i = 0; i < 3; i++)
		{
			printf("%s: ", labels[i]);
			for (j = 0; j < ITERATIONS; j++)
			{
				if (j != 0) printf(", ");

				printf("%ld", (long int)((double)(CPG)*(double)(GIG * time_stamp[i][j].tv_sec + time_stamp[i][j].tv_nsec)));
			}
		printf("\n");
		}


		delete iimg;
		delete ipts;

	}

	else
	{
		printf("Invalid use!\n");
	}

	return 0;
}

/*************************************************/

void PaintSURF(bitmap_image &img, std::vector<IPoint> * ipts)
{
	image_drawer idrawer(img);

	for (std::vector<IPoint>::iterator ip = ipts->begin(); ip != ipts->end(); ++ip)
	{
		unsigned char myPen[3] = {0x00, 0x00, 0x00};
		int S = 2 * (int)((float)2.5 * ip->scale);
		int R = (int)(S / (float)2.0);

		int  pt[2] = {(int)(ip->x), (int)(ip->y)};
		int ptR[2] = {(int)(R * cos(ip->orientation)),
			      (int)(R * sin(ip->orientation))};

		myPen[(ip->laplacian > 0 ? 2 : 0)] = 0xFF;
		idrawer.pen_color(myPen[0], myPen[1], myPen[2]);
		idrawer.circle(pt[0] - R, pt[1] - R, R);

		idrawer.pen_color(0x00, 0xFF, 0x00);
		idrawer.line_segment(pt[0], pt[1], pt[0] + ptR[0], pt[1] + ptR[1]);
      }
}