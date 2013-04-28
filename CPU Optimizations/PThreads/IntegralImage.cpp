#include "IntegralImage.h"
#include "bitmap_image.hpp"
#include <pthread.h>
#include <cmath>
#include "parameters.h"

pthread_barrier_t barr_ii;
struct thread_data {
	int thread_id;
	IntegralImage * pic;
	bitmap_image * image;
};


const float IntegralImage::cR = (float).2989;
const float IntegralImage::cG = (float).5870;
const float IntegralImage::cB = (float).1140;


float IntegralImage::getValue(int y, int x)
{
//	return this->Matrix[y*Height+x];
	return this->Matrix[y][x];
}

void IntegralImage::setValue(int y, int x, float value)
{ 
//	this->Matrix[y*Height+x] = value;
	this->Matrix[y][x] = value;
}


IntegralImage::IntegralImage()
{
	this->Width = 0;
	this->Height = 0;
	this->Matrix = NULL;
}

IntegralImage::IntegralImage(int width, int height)
{
	this->Width = width;
	this->Height = height;

//	this->Matrix = (float*)malloc(sizeof(float)*width*height);

	this->Matrix = new float*[height];
	for (int ii = 0; ii < height; ii++)
		Matrix[ii] = new float[width];

}

IntegralImage::~IntegralImage()
{
	if (this->Matrix != NULL) {
		for (int ii = 0; ii < Height; ii++)
			delete [] Matrix[ii];
		delete [] Matrix;
	}
	
//	delete [] Matrix;

}

IntegralImage * IntegralImage::FromImage(bitmap_image &image)
{
	void  * FromImage_work(void * threadarg);
	IntegralImage * pic = new IntegralImage(image.width(), image.height());
	pthread_t threads[NUM_THREADS];
	struct thread_data data[NUM_THREADS];

	if (pthread_barrier_init(&barr_ii, NULL, NUM_THREADS))
	{
		printf("Could not createa barrier\n");
		exit(-1);
	}

	for (int t = 0; t < NUM_THREADS; t++)
	{
		data[t].thread_id = t;
		data[t].pic = pic;
		data[t].image = &image;

		int rc = pthread_create(&threads[t], NULL, FromImage_work, (void*) &(data[t]));
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
			printf("\n ERROR on join\n");
			exit(19);
		}
	}

	return pic;
}

void  * FromImage_work(void * threadarg)
{
	struct thread_data * my_data = (struct thread_data *) threadarg;
	IntegralImage * pic = my_data->pic;
	bitmap_image * image = my_data->image;
	int thread_id = my_data-> thread_id;


	float colsum = 0;
	float rowsum = 0;
	unsigned char red, green, blue;

	unsigned int start = (image->height() / NUM_THREADS) * thread_id;
	unsigned int end = (image->height() / NUM_THREADS) + start;
	if (end > image->height()) end  = image->height();

	for (unsigned int y = start; y < end; y++)
	{
		rowsum = 0;
		for (unsigned int x = 0; x < image->width(); x++)
		{
			image->get_pixel(x, y, red, green, blue);
			rowsum += (pic->cR * red + pic->cG * green + pic->cB * blue) / (float)255;

			pic->setValue(y, x, rowsum);
		}
	}


	int rc = pthread_barrier_wait(&barr_ii);
	if (rc != 0 && rc!= PTHREAD_BARRIER_SERIAL_THREAD)
	{
		printf("Could not wait on barrier\n");
		exit(-1);
	}

	start = (image->width() / NUM_THREADS) * thread_id;
	end = (image->width() / NUM_THREADS) + start;
	if (end > image->width()) end  = image->width();

	for (unsigned int x = start; x < end; x++)
	{
		colsum = pic->getValue(0, x);
		for (unsigned int y = 1; y < image->height(); y++)
		{
			colsum += pic->getValue(y,x);
			pic->setValue(y, x, colsum);
		}
	}

	pthread_exit(NULL);
}


float IntegralImage::BoxIntegral(int row, int col, int rows, int cols)
{
	//printf("row = %d, Height = %d, col = %d, Width = %d\n", row,Height,col,Width);
	// The subtraction by one for row/col is because row/col is inclusive.
	int r1 = std::min(row, Height) - 1;
	int c1 = std::min(col, Width) - 1;
	int r2 = std::min(row + rows, Height) - 1;
	int c2 = std::min(col + cols, Width) - 1;
	//printf("After min\n");
	float A = 0, B = 0, C = 0, D = 0;

	//printf("r1 = %d, c1 = %d, r2 = %d, c2 =%d", r1,c1,r2,r2);
/*
	if (r1 >= 0 && c1 >= 0) A = Matrix[r1*Height+c1];
	if (r1 >= 0 && c2 >= 0) B = Matrix[r1*Height+c2];
	if (r2 >= 0 && c1 >= 0) C = Matrix[r2*Height+c1];
	if (r2 >= 0 && c2 >= 0) D = Matrix[r2*Height+c2];
*/
	if (r1 >= 0 && c1 >= 0) A = Matrix[r1][c1];
	if (r1 >= 0 && c2 >= 0) B = Matrix[r1][c2];
	if (r2 >= 0 && c1 >= 0) C = Matrix[r2][c1];
	if (r2 >= 0 && c2 >= 0) D = Matrix[r2][c2];

	return std::max((float)0, A - B - C + D);
}

/// <summary>
/// Get Haar Wavelet X repsonse
/// </summary>
/// <param name="row"></param>
/// <param name="column"></param>
/// <param name="size"></param>
/// <returns></returns>
float IntegralImage::HaarX(int row, int column, int size)
{
	return BoxIntegral(row - size / 2, column, size, size / 2)
	- 1 * BoxIntegral(row - size / 2, column - size / 2, size, size / 2);
}

/// <summary>
/// Get Haar Wavelet Y repsonse
/// </summary>
/// <param name="row"></param>
/// <param name="column"></param>
/// <param name="size"></param>
/// <returns></returns>
float IntegralImage::HaarY(int row, int column, int size)
{
	return BoxIntegral(row, column - size / 2, size / 2, size)
	- 1 * BoxIntegral(row - size / 2, column - size / 2, size / 2, size);
}

