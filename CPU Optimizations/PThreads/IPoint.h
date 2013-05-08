#ifndef IPOINT_H
#define IPOINT_H

class IPoint
{

public:

	IPoint()
	{
		orientation = 0;
		descriptor = 0;
	}
	
	~IPoint()
	{
		if (descriptor != 0)
			delete [] descriptor;
	}

	// Coordinates of the detected interest point
	float x, y;

	// Detected scale
	float scale;

	// Response of the detected feature (strength)
	float response;

	// Orientation measured anti-clockwise from +ve x-axis
	float orientation;

	// Sign of laplacian for fast matching purposes
	int laplacian;

	// Descriptor vector
	int descriptorLength;
	float * descriptor;

	void SetDescriptorLength(int Size)
	{
		descriptorLength = Size;
		if (descriptor != 0)
			delete [] descriptor;
		descriptor = new float[Size];
	}
};

#endif
