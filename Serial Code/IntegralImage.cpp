#include <algorithm>
#include "IntegralImage.h"
#include "bitmap_image.hpp"


float getValue(int y, int x)
{
    return this.Matrix[y][x]
}

void setValue(int y, int x)
{ 
	this.Matrix[y][x] = value; 
}

IntegralImage(int width, int height)
{
    this.Width = width;
    this.Height = height;

    this.Matrix = new float* [height];
	for(int i = 0; i < height; i++)
	{
		this.Matrix[i] = new float* [width];
	}
}

public static IntegralImage FromImage(bitmap_image image)
{
    IntegralImage pic = new IntegralImage(image.width, image.height);

    float rowsum = 0;
	unsigned char red, green, blue;

    for (int x = 0; x < image.width; x++)
    {
		image.GetPixel(x, 0, &red, &green, &blue);
		rowsum += (cR * red + cG * green + cB * blue) / 255f;
		pic[0][x] = rowsum;
    }


    for (int y = 1; y < image.height; y++)
    {
		rowsum = 0;
		for (int x = 0; x < image.Width; x++)
		{
			image.GetPixel(x, y, &red, &green, &blue);
			rowsum += (cR * red + cG * green + cB * blue) / 255f;

			// integral image is rowsum + value above        
			pic[y][x] = rowsum + pic[y - 1][x];
		}
    }

    return pic;
}


public float BoxIntegral(int row, int col, int rows, int cols)
{
    // The subtraction by one for row/col is because row/col is inclusive.
    int r1 = min(row, Height) - 1;
    int c1 = min(col, Width) - 1;
    int r2 = min(row + rows, Height) - 1;
    int c2 = min(col + cols, Width) - 1;

    float A = 0, B = 0, C = 0, D = 0;
    if (r1 >= 0 && c1 >= 0) A = Matrix[r1][c1];
    if (r1 >= 0 && c2 >= 0) B = Matrix[r1][c2];
    if (r2 >= 0 && c1 >= 0) C = Matrix[r2][c1];
    if (r2 >= 0 && c2 >= 0) D = Matrix[r2][c2];

    return max(0, A - B - C + D);
}

/// <summary>
/// Get Haar Wavelet X repsonse
/// </summary>
/// <param name="row"></param>
/// <param name="column"></param>
/// <param name="size"></param>
/// <returns></returns>
public float HaarX(int row, int column, int size)
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
public float HaarY(int row, int column, int size)
{
    return BoxIntegral(row, column - size / 2, size / 2, size)
    - 1 * BoxIntegral(row - size / 2, column - size / 2, size / 2, size);
}
}
