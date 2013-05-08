#include <iostream>
#include "FastHessian.h"
#include "IntegralImage.h"
#include "bitmap_image.hpp"
#include "IPoint.h"
#include "SurfDescriptor.h"
#include <vector>
#include <cmath>
#include <time.h>

#define CPG 2.53           // Cycles per GHz -- Adjust to your computer
#define ITERATIONS 10
#define TOL 1

#define GIG 1000000000

/************************************************************************************
*  Usage:
*      ./surf [source file name] [destination file name]
*
************************************************************************************/

struct timespec diff(struct timespec start, struct timespec end);

int main(int argc, char ** argv)
{
	if (argc == 3)
	{

		int tests = 0;
		void PaintSURF(bitmap_image &img, std::vector<IPoint> * ipts);
		int clock_gettime(clockid_t clk_id, struct timespec *tp);

		struct timespec time1, time2;
		struct timespec time_stamp_all[ITERATIONS];
		struct timespec time_stamp[3][ITERATIONS];

		long int i, j;

		std::string file_name(argv[1]);
		bitmap_image image(file_name);
		printf("File: %s, Width: %d, Height: %d\n", argv[1], image.width(), image.height());
		IntegralImage * test_img = IntegralImage::FromImage(image);
		delete test_img;

		IntegralImage * iimg;
		std::vector<IPoint> * ipts = NULL;

		for (int ii = 0; ii < ITERATIONS; ii++)
		{
			clock_gettime(CLOCK_REALTIME, &time1);

			iimg = IntegralImage::FromImage(image);

			ipts = NULL;
			ipts = FastHessian::getIpoints((float)0.0002, 5, 2, iimg);

			if (ipts == NULL)
			{
				printf("getIpoints failed\n");
				return 0;
			}

			SurfDescriptor::DecribeInterestPoints(ipts, false, false, iimg);
			clock_gettime(CLOCK_REALTIME, &time2);
			time_stamp_all[ii] = diff(time1,time2);

			if (ii != ITERATIONS - 1)
			{
				delete iimg;
				delete ipts;
			}
		}

		printf("Detected Features: %d\n", ipts->size());

		printf("Cumulative Time: ");
		for (j = 0; j < ITERATIONS; j++)
		{
			if (j != 0) printf(", ");

			printf("%.2f", ((double)(GIG * time_stamp_all[j].tv_sec + time_stamp_all[j].tv_nsec) / (double)(1000000)));
		}
		printf("\n");

		delete iimg;
		delete ipts;

		const char * labels[3] = {"IntegralImage","FastHessian","SurfDescriptor"};

		for (int ii = 0; ii < ITERATIONS; ii++)
		{
			clock_gettime(CLOCK_REALTIME, &time1);
			iimg = IntegralImage::FromImage(image);
			clock_gettime(CLOCK_REALTIME, &time2);
			time_stamp[tests][ii] = diff(time1,time2);
			if (ii != ITERATIONS - 1) delete iimg;
		}
		tests++;

		ipts = NULL;

		for (int ii = 0; ii < ITERATIONS; ii++)
		{
			clock_gettime(CLOCK_REALTIME, &time1);
			ipts = FastHessian::getIpoints((float)0.0002, 5, 2, iimg);
			clock_gettime(CLOCK_REALTIME, &time2);
			time_stamp[tests][ii] = diff(time1,time2);
			if (ii != ITERATIONS - 1) delete ipts;
		}
		tests++;


		if (ipts == NULL)
		{
			printf("getIpoints failed\n");
			return 0;
		}

		for (int ii = 0; ii < ITERATIONS; ii++)
		{
			clock_gettime(CLOCK_REALTIME, &time1);
			SurfDescriptor::DecribeInterestPoints(ipts, false, false, iimg);
			clock_gettime(CLOCK_REALTIME, &time2);
			time_stamp[tests][ii] = diff(time1,time2);
		}
		tests++;

		PaintSURF(image, ipts);
		image.save_image(argv[2]);

		for (i = 0; i < tests; i++)
		{
			printf("%s: ", labels[i]);
			for (j = 0; j < ITERATIONS; j++)
			{
				if (j != 0) printf(", ");

				printf("%.2f", ((double)(GIG * time_stamp[i][j].tv_sec + time_stamp[i][j].tv_nsec) / (double)(1000000)));
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

struct timespec diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

/*************************************************/

void PaintSURF(bitmap_image &img, std::vector<IPoint> * ipts)
{

	image_drawer idrawer(img);

	int count;
	for (std::vector<IPoint>::iterator ip = ipts->begin(); ip != ipts->end(); ++ip, count++)
	{
//		printf("Point %d: (%f, %f), scale: %f, response %f, orientation: %f, laplacian: %d\n", count, ip->x, ip->y, ip->scale, ip->response, ip->orientation, ip->laplacian);
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
