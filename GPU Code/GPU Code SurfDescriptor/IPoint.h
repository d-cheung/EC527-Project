#ifndef IPOINT_H
#define IPOINT_H

class IPoint
{

public:

	/// <summary>
	/// Default ctor
	/// </summary>
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

	/// <summary>
	/// Coordinates of the detected interest point
	/// </summary>
	float x, y;

	/// <summary>
	/// Detected scale
	/// </summary>
	float scale;

	/// <summary>
	/// Response of the detected feature (strength)
	/// </summary>
	float response;

	/// <summary>
	/// Orientation measured anti-clockwise from +ve x-axis
	/// </summary>
	float orientation;

	/// <summary>
	/// Sign of laplacian for fast matching purposes
	/// </summary>
	int laplacian;

	/// <summary>
	/// Descriptor vector
	/// </summary>
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
