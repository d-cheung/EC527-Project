#include "IPoint.h"
#include <cmath>
#include <vector>
#include "FastHessian.h"
#include "IntegralImage.h"

// Static one-call do it all method
std::vector<IPoint> * FastHessian::getIpoints(float thresh, int octaves, int init_sample, IntegralImage * img)
{
	FastHessian fh(thresh, octaves, init_sample, img);
	return fh.getIpoints();
}


// Constructor with parameters
FastHessian::FastHessian(float thresh, int octaves, int init_sample, IntegralImage * img)
{
	this->thresh = thresh;
	this->octaves = octaves;
	this->init_sample = init_sample;
	this->img = img;
	this->ipts = NULL;
	this->responseMap = NULL;
}

FastHessian::~FastHessian()
{
	if (responseMap != NULL) {
		for (int ii = 0; ii < responseMap->size(); ii++)
			delete (*responseMap)[ii];
		delete responseMap;
	}
}

// Find the image features and write into vector of features
std::vector<IPoint> * FastHessian::getIpoints()
{
	// filter index map
	int filter_map[5][4] = {{0,1,2,3}, {1,3,4,5}, {3,5,6,7}, {5,7,8,9}, {7,9,10,11}};

	// Clear the vector of exisiting ipts
	if (ipts == NULL) ipts = new std::vector<IPoint>();
	else ipts->clear();

	// Build the response map
	buildResponseMap();

	// Get the response layers
	ResponseLayer * b, * m, * t;
	for (int o = 0; o < octaves; ++o) for (int i = 0; i <= 1; ++i)
	{
		b = (*responseMap)[filter_map[o][i]];
		m = (*responseMap)[filter_map[o][i+1]];
		t = (*responseMap)[filter_map[o][i+2]];

		// loop over middle response layer at density of the most 
		// sparse layer (always top), to find maxima across scale and space
		for (int r = 0; r < t->height; ++r)
		{
			for (int c = 0; c < t->width; ++c)
			{
				if (isExtremum(r, c, *t, *m, *b))
				{
					interpolateExtremum(r, c, *t, *m, *b);
				}
			}
		}
	}
	return ipts;
}


// Build map of DoH responses
void FastHessian::buildResponseMap()
{
	// Calculate responses for the first 4 octaves:
	// Oct1: 9,  15, 21, 27
	// Oct2: 15, 27, 39, 51
	// Oct3: 27, 51, 75, 99
	// Oct4: 51, 99, 147,195
	// Oct5: 99, 195,291,387

	// Deallocate memory and clear any existing response layers
	if (responseMap == NULL) responseMap = new std::vector<ResponseLayer *>();
	else responseMap->clear();

	// Get image attributes
	int w = (img->Width / init_sample);
	int h = (img->Height / init_sample);
	int s = (init_sample);

	// Calculate approximated determinant of hessian values
	if (octaves >= 1)
	{
		responseMap->push_back(new ResponseLayer(w, h, s,  9));
		responseMap->push_back(new ResponseLayer(w, h, s, 15));
		responseMap->push_back(new ResponseLayer(w, h, s, 21));
		responseMap->push_back(new ResponseLayer(w, h, s, 27));
	}
	 
	if (octaves >= 2)
	{
		responseMap->push_back(new ResponseLayer(w / 2, h / 2, s * 2, 39));
		responseMap->push_back(new ResponseLayer(w / 2, h / 2, s * 2, 51));
	}

	if (octaves >= 3)
	{
		responseMap->push_back(new ResponseLayer(w / 4, h / 4, s * 4, 75));
		responseMap->push_back(new ResponseLayer(w / 4, h / 4, s * 4, 99));
	}

	if (octaves >= 4)
	{
		responseMap->push_back(new ResponseLayer(w / 8, h / 8, s * 8, 147));
		responseMap->push_back(new ResponseLayer(w / 8, h / 8, s * 8, 195));
	}

	if (octaves >= 5)
	{
		responseMap->push_back(new ResponseLayer(w / 16, h / 16, s * 16, 291));
		responseMap->push_back(new ResponseLayer(w / 16, h / 16, s * 16, 387));
	}

	  // Extract responses from the image
	for (unsigned int i = 0; i < responseMap->size(); ++i)
	{
		buildResponseLayer(*((*responseMap)[i]));
	}
}


// Build Responses for a given ResponseLayer
void FastHessian::buildResponseLayer(ResponseLayer &rl)
{
	int step = rl.step;                      // step size for this filter
	int b = (rl.filter - 1) / 2;             // border for this filter
	int l = rl.filter / 3;                   // lobe for this filter (filter size / 3)
	int w = rl.filter;                       // filter size
	float inverse_area = (float)1.0 / (w * w);       // normalisation factor
	float Dxx0, Dyy0, Dxy0;
	float Dxx1, Dyy1, Dxy1;
	float Dxx2, Dyy2, Dxy2;
	float Dxx3, Dyy3, Dxy3;

	for (int r, c0, c1, c2 ,c3, ar = 0, index = 0; ar < rl.height; ++ar)
	{
		int ac;
		for (ac = 0; ac < rl.width-4; ac+=4, index+=4)
		{
			// get the image coordinates
			r = ar * step;
			c0 = (ac+0) * step;
			c1 = (ac+1) * step;
			c2 = (ac+2) * step;
			c3 = (ac+3) * step;

			// Compute response components
			Dxx0 = img->BoxIntegral(r - l + 1, (c0) - b, 2 * l - 1, w)
			     - img->BoxIntegral(r - l + 1, (c0) - l / 2, 2 * l - 1, l) * 3;
			Dxx1 = img->BoxIntegral(r - l + 1, (c1) - b, 2 * l - 1, w)
			     - img->BoxIntegral(r - l + 1, (c1) - l / 2, 2 * l - 1, l) * 3;
			Dxx2 = img->BoxIntegral(r - l + 1, (c2) - b, 2 * l - 1, w)
			     - img->BoxIntegral(r - l + 1, (c2) - l / 2, 2 * l - 1, l) * 3;
			Dxx3 = img->BoxIntegral(r - l + 1, (c3) - b, 2 * l - 1, w)
			     - img->BoxIntegral(r - l + 1, (c3) - l / 2, 2 * l - 1, l) * 3;

			Dyy0 = img->BoxIntegral(r - b, (c0) - l + 1, w, 2 * l - 1)
			     - img->BoxIntegral(r - l / 2, (c0) - l + 1, l, 2 * l - 1) * 3;
			Dyy1 = img->BoxIntegral(r - b, (c1) - l + 1, w, 2 * l - 1)
			     - img->BoxIntegral(r - l / 2, (c1) - l + 1, l, 2 * l - 1) * 3;
			Dyy2 = img->BoxIntegral(r - b, (c2) - l + 1, w, 2 * l - 1)
			     - img->BoxIntegral(r - l / 2, (c2) - l + 1, l, 2 * l - 1) * 3;
			Dyy3 = img->BoxIntegral(r - b, (c3) - l + 1, w, 2 * l - 1)
			     - img->BoxIntegral(r - l / 2, (c3) - l + 1, l, 2 * l - 1) * 3;

			Dxy0 = + img->BoxIntegral(r - l, (c0) + 1, l, l)
			       + img->BoxIntegral(r + 1, (c0) - l, l, l)
			       - img->BoxIntegral(r - l, (c0) - l, l, l)
			       - img->BoxIntegral(r + 1, (c0) + 1, l, l);
			Dxy1 = + img->BoxIntegral(r - l, (c0) + 1, l, l)
			       + img->BoxIntegral(r + 1, (c1) - l, l, l)
			       - img->BoxIntegral(r - l, (c2) - l, l, l)
			       - img->BoxIntegral(r + 1, (c3) + 1, l, l);
			Dxy2 = + img->BoxIntegral(r - l, (c0) + 1, l, l)
			       + img->BoxIntegral(r + 1, (c1) - l, l, l)
			       - img->BoxIntegral(r - l, (c2) - l, l, l)
			       - img->BoxIntegral(r + 1, (c3) + 1, l, l);
			Dxy3 = + img->BoxIntegral(r - l, (c0) + 1, l, l)
			       + img->BoxIntegral(r + 1, (c1) - l, l, l)
			       - img->BoxIntegral(r - l, (c2) - l, l, l)
			       - img->BoxIntegral(r + 1, (c3) + 1, l, l);

			// Normalise the filter responses with respect to their size
			Dxx0 *= inverse_area;
			Dxx1 *= inverse_area;
			Dxx2 *= inverse_area;
			Dxx3 *= inverse_area;

			Dyy0 *= inverse_area;
			Dyy1 *= inverse_area;
			Dyy2 *= inverse_area;
			Dyy3 *= inverse_area;

			Dxy0 *= inverse_area;
			Dxy1 *= inverse_area;
			Dxy2 *= inverse_area;
			Dxy3 *= inverse_area;

			// Get the determinant of hessian response & laplacian sign
			rl.responses[index+0] = (Dxx0 * Dyy0 - (float)0.81 * Dxy0 * Dxy0);
			rl.responses[index+1] = (Dxx1 * Dyy1 - (float)0.81 * Dxy1 * Dxy1);
			rl.responses[index+2] = (Dxx2 * Dyy2 - (float)0.81 * Dxy2 * Dxy2);
			rl.responses[index+3] = (Dxx3 * Dyy3 - (float)0.81 * Dxy3 * Dxy3);

			rl.laplacian[index+0] = (unsigned char)(Dxx0 + Dyy0 >= 0 ? 1 : 0);
			rl.laplacian[index+1] = (unsigned char)(Dxx1 + Dyy1 >= 0 ? 1 : 0);
			rl.laplacian[index+2] = (unsigned char)(Dxx2 + Dyy2 >= 0 ? 1 : 0);
			rl.laplacian[index+3] = (unsigned char)(Dxx3 + Dyy3 >= 0 ? 1 : 0);
		}

		for (; ac < rl.width; ++ac, index++)
		{
			// get the image coordinates
			r = ar * step;
			c0 = (ac+0) * step;

			// Compute response components
			Dxx0 = img->BoxIntegral(r - l + 1, (c0) - b, 2 * l - 1, w)
			     - img->BoxIntegral(r - l + 1, (c0) - l / 2, 2 * l - 1, l) * 3;
			Dyy0 = img->BoxIntegral(r - b, (c0) - l + 1, w, 2 * l - 1)
			     - img->BoxIntegral(r - l / 2, (c0) - l + 1, l, 2 * l - 1) * 3;
			Dxy0 = + img->BoxIntegral(r - l, (c0) + 1, l, l)
			       + img->BoxIntegral(r + 1, (c0) - l, l, l)
			       - img->BoxIntegral(r - l, (c0) - l, l, l)
			       - img->BoxIntegral(r + 1, (c0) + 1, l, l);

			// Normalise the filter responses with respect to their size
			Dxx0 *= inverse_area;
			Dyy0 *= inverse_area;
			Dxy0 *= inverse_area;

			// Get the determinant of hessian response & laplacian sign
			rl.responses[index+0] = (Dxx0 * Dyy0 - (float)0.81 * Dxy0 * Dxy0);
			rl.laplacian[index+0] = (unsigned char)(Dxx0 + Dyy0 >= 0 ? 1 : 0);
		}
	}
}


// Test whether the point r,c in the middle layer is extremum in 3x3x3 neighbourhood
bool FastHessian::isExtremum(int r, int c, ResponseLayer &t, ResponseLayer &m, ResponseLayer &b)
{
	// bounds check
	int layerBorder = (t.filter + 1) / (2 * t.step);
	if (r <= layerBorder || r >= t.height - layerBorder || c <= layerBorder || c >= t.width - layerBorder)
		return false;

	// check the candidate point in the middle layer is above thresh 
	float candidate = m.getResponse(r, c, t);
	if (candidate < thresh)
		return false;

	for (int rr = -1; rr <= 1; ++rr)
	{
		for (int cc = -1; cc <= 1; ++cc)
		{
			// if any response in 3x3x3 is greater candidate not maximum
			if (t.getResponse(r + rr, c + cc) >= candidate ||
			    ((rr != 0 || cc != 0) && m.getResponse(r + rr, c + cc, t) >= candidate) ||
				b.getResponse(r + rr, c + cc, t) >= candidate)
			{
				return false;
			}
		}
	}

	return true;
}


// Interpolate scale-space extrema to subpixel accuracy to form an image feature
void FastHessian::interpolateExtremum(int r, int c, ResponseLayer &t, ResponseLayer &m, ResponseLayer &b)
{
	double * D  = BuildDerivative(r, c, t, m, b);
	double * H  = BuildHessian(r, c, t, m, b);
	double * Hi = Inverse(H);

	if (Hi != NULL) {
		double * Of = MMM_neg_3x3_3x1(Hi, D);

		// get the offsets from the interpolation
		double O[3] = { Of[0], Of[1], Of[2] };

		// get the step distance between filters
		int filterStep = (m.filter - b.filter);
	 
		// If point is sufficiently close to the actual extremum
		if (fabs(O[0]) < (float)0.5 && fabs(O[1]) < (float)0.5 && fabs(O[2]) < (float)0.5)
		{
			IPoint ipt;
			ipt.x = (float)((c + O[0]) * t.step);
			ipt.y = (float)((r + O[1]) * t.step);
			ipt.scale = (float)((0.1333f) * (m.filter + O[2] * filterStep));
			ipt.laplacian = (int)(m.getLaplacian(r,c,t));
			ipts->push_back(ipt);
		}

		delete[] Of;
		delete[] Hi;
	}


	delete[] D;
	delete[] H;

}

// Build Matrix of First Order Scale-Space derivatives
double * FastHessian::BuildDerivative(int r, int c, ResponseLayer &t, ResponseLayer &m, ResponseLayer &b)
{
	double dx, dy, ds;

	dx = (m.getResponse(r, c + 1, t) - m.getResponse(r, c - 1, t)) / (float)(2.0);
	dy = (m.getResponse(r + 1, c, t) - m.getResponse(r - 1, c, t)) / (float)(2.0);
	ds = (t.getResponse(r, c)        - b.getResponse(r, c, t))     / (float)(2.0);

	double * D = new double[3];
	D[0] = dx;
	D[1] = dy;
	D[2] = ds;
	return D;
}


// Build Hessian Matrix 
double * FastHessian::BuildHessian(int r, int c, ResponseLayer &t, ResponseLayer &m, ResponseLayer &b)
{
	double v, dxx, dyy, dss, dxy, dxs, dys;

	v = m.getResponse(r, c, t);
	dxx = m.getResponse(r, c + 1, t) + m.getResponse(r, c - 1, t) - 2 * v;
	dyy = m.getResponse(r + 1, c, t) + m.getResponse(r - 1, c, t) - 2 * v;
	dss = t.getResponse(r, c) + b.getResponse(r, c, t) - 2 * v;
	dxy = (m.getResponse(r + 1, c + 1, t) - m.getResponse(r + 1, c - 1, t) -
	       m.getResponse(r - 1, c + 1, t) + m.getResponse(r - 1, c - 1, t)) / (float)(4.0);
	dxs = (t.getResponse(r, c + 1) - t.getResponse(r, c - 1) -
	       b.getResponse(r, c + 1, t) + b.getResponse(r, c - 1, t)) / (float)(4.0);
	dys = (t.getResponse(r + 1, c) - t.getResponse(r - 1, c) -
	       b.getResponse(r + 1, c, t) + b.getResponse(r - 1, c, t)) / (float)(4.0);

	double * H = new double[9];
	H[0] = dxx;
	H[1] = dxy;
	H[2] = dxs;
	H[3] = dxy;
	H[4] = dyy;
	H[5] = dys;
	H[6] = dxs;
	H[7] = dys;
	H[8] = dss;
	return H;
}

// Return inverse of the Matrix m
double * FastHessian::Inverse(double * m)
{
	double det = ((m[0]*m[4]*m[8]) + (m[1]*m[5]*m[6]) + (m[2]*m[3]*m[7]) -
		      (m[2]*m[4]*m[6]) - (m[1]*m[3]*m[8]) - (m[0]*m[5]*m[7]));

	if (det == 0) return NULL;

	double * invert = new double[9];
	double A =      (m[4] * m[8] - m[5] * m[7]) / det;
	double B = -1 * (m[3] * m[8] - m[6] * m[5]) / det;
	double C =      (m[3] * m[7] - m[4] * m[6]) / det;
	double D = -1 * (m[1] * m[8] - m[2] * m[7]) / det;
	double E =      (m[0] * m[8] - m[2] * m[6]) / det;
	double F = -1 * (m[0] * m[7] - m[1] * m[6]) / det;
	double G =      (m[1] * m[5] - m[2] * m[4]) / det;
	double H = -1 * (m[0] * m[5] - m[2] * m[3]) / det;
	double K =      (m[0] * m[4] - m[1] * m[3]) / det;
	
	invert[0] = A;
	invert[1] = D;
	invert[2] = G;
	invert[3] = B;
	invert[4] = E;
	invert[5] = H;
	invert[6] = C;
	invert[7] = F;
	invert[8] = K;

	return invert;

}

// Perform -1 * (A * B) where A is a 3x3 matrix and B is a 3x1 matrix
double * FastHessian::MMM_neg_3x3_3x1(double * A, double * B)
{
	double * C = new double[3];
	double a = -1 * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
	double b = -1 * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
	double c = -1 * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);

	C[0] = a;
	C[1] = b;
	C[2] = c;
	
	return C;
}


FastHessian::ResponseLayer::ResponseLayer()
{
	this->width = 0;
	this->height = 0;
	this->step = 0;
	this->filter = 0;
	this->responses = NULL;
	this->laplacian = NULL;
}

FastHessian::ResponseLayer::ResponseLayer(int width, int height, int step, int filter)
{
	this->width = width;
	this->height = height;
	this->step = step;
	this->filter = filter;

	responses = new float[width * height];
	laplacian = new unsigned char[width * height];
}

FastHessian::ResponseLayer::~ResponseLayer()
{
	if (responses != NULL)
		delete [] responses;
	if (laplacian != NULL)
		delete [] laplacian;
}

unsigned char FastHessian::ResponseLayer::getLaplacian(int row, int column)
{
	return laplacian[row * width + column];
}

unsigned char FastHessian::ResponseLayer::getLaplacian(int row, int column, ResponseLayer &src)
{
	int scale = this->width / src.width;
	return laplacian[(scale * row) * width + (scale * column)];
}

float FastHessian::ResponseLayer::getResponse(int row, int column)
{
	return responses[row * width + column];
}

float FastHessian::ResponseLayer::getResponse(int row, int column, ResponseLayer &src)
{
	int scale = this->width / src.width;
	return responses[(scale * row) * width + (scale * column)];
}
