#include <iostream>
#include "FastHessian.h"
#include "IntegralImage.h"
#include "bitmap_image.hpp"
#include "IPoint.h"
#include "SurfDescriptor.h"
#include <vector>
#include <cmath>

using namespace std;

void PaintSURF(bitmap_image &img, vector<IPoint> * ipts);

int main(int argc, char ** argv)
{
	if (argc == 3)
	{

		string file_name(argv[1]);
		bitmap_image image(file_name);
//		image.save_image("test01_saved.bmp");

		printf("%d x %d\n",image.width(), image.height());
		IntegralImage * iimg = IntegralImage::FromImage(image);

		vector<IPoint> * ipts = NULL;
		ipts = FastHessian::getIpoints((float)0.0002, 5, 2, iimg);

		if (ipts == NULL)
		{
			printf("getIpoints failed\n");
			return 0;
		}

		SurfDescriptor::DecribeInterestPoints(ipts, false, false, iimg);
		PaintSURF(image, ipts);
		image.save_image(argv[2]);

		delete iimg;
		delete ipts;

	}

	else
	{
		printf("Invalid use!  Need bitmap file\n");
	}

	return 0;
}


void PaintSURF(bitmap_image &img, vector<IPoint> * ipts)
{

	unsigned char   red[3] = {0xFF, 0x00, 0x00};
	unsigned char green[3] = {0x00, 0xFF, 0x00};
	unsigned char  blue[3] = {0x00, 0x00, 0xFF};
	image_drawer idrawer(img);

	for (std::vector<IPoint>::iterator ip = ipts->begin(); ip != ipts->end(); ++ip)
	{
		unsigned char myPen[3] = {0x00, 0x00, 0x00};
		int S = 2 * (int)((float)2.5 * ip->scale);
		int R = (int)(S / (float)2.0);

		int  pt[2] = {(int)(ip->x), (int)(ip->y)};
		int ptR[2] = {(int)(R * cos(ip->orientation)), (int)(R * sin(ip->orientation))};

		myPen[(ip->laplacian > 0 ? 2 : 0)] = 0xFF;
		idrawer.pen_color(myPen[0], myPen[1], myPen[2]);
		idrawer.circle(pt[0] - R, pt[1] - R, R);

		idrawer.pen_color(0x00, 0xFF, 0x00);
		idrawer.line_segment(pt[0], pt[1], pt[0] + ptR[0], pt[1] + ptR[1]);
      }
}
