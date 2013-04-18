class IPoint
{

public:

    /// <summary>
    /// Default ctor
    /// </summary>
    IPoint()
    {
		orientation = 0;
    }
	
	~IPoint()
	{
		delete[] descriptor;
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
    float [] descriptor = NULL;
    void SetDescriptorLength(int Size)
    {
		descriptorLength = Size;
		descriptor = new float[Size];
    }
  }
}
