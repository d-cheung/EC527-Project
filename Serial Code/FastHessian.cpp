#include "IPoint.h"
#include <cmath>
#include <vector>
#include "FastHessian.h"
#include "IntegralImage.h"

/// <summary>
/// Static one-call do it all method
/// </summary>
/// <param name="thresh"></param>
/// <param name="octaves"></param>
/// <param name="init_sample"></param>
/// <param name="img"></param>
/// <returns></returns>
std::vector<IPoint> * FastHessian::getIpoints(float thresh, int octaves, int init_sample, IntegralImage img)
{
	FastHessian * fh = new FastHessian(thresh, octaves, init_sample, img);
	return fh->getIpoints();
}


/// <summary>
/// Constructor with parameters
/// </summary>
/// <param name="thresh"></param>
/// <param name="octaves"></param>
/// <param name="init_sample"></param>
/// <param name="img"></param>
FastHessian::FastHessian(float thresh, int octaves, int init_sample, IntegralImage img)
{
	this->thresh = thresh;
	this->octaves = octaves;
	this->init_sample = init_sample;
	this->img = img;
}

/// <summary>
/// Find the image features and write into vector of features
/// </summary>
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
	ResponseLayer b, m, t;
	for (int o = 0; o < octaves; ++o) for (int i = 0; i <= 1; ++i)
	{
		b = (*responseMap)[filter_map[o][i]];
		m = (*responseMap)[filter_map[o][i+1]];
		t = (*responseMap)[filter_map[o][i+2]];

		// loop over middle response layer at density of the most 
		// sparse layer (always top), to find maxima across scale and space
		for (int r = 0; r < t.height; ++r)
		{
			for (int c = 0; c < t.width; ++c)
			{
				if (isExtremum(r, c, t, m, b))
				{
					interpolateExtremum(r, c, t, m, b);
				}
			}
		}
	}

	return ipts;
}


/// <summary>
/// Build map of DoH responses
/// </summary>
void FastHessian::buildResponseMap()
{
	// Calculate responses for the first 4 octaves:
	// Oct1: 9,  15, 21, 27
	// Oct2: 15, 27, 39, 51
	// Oct3: 27, 51, 75, 99
	// Oct4: 51, 99, 147,195
	// Oct5: 99, 195,291,387

	// Deallocate memory and clear any existing response layers
	if (responseMap == NULL) responseMap = new std::vector<ResponseLayer>();
	else responseMap->clear();

	// Get image attributes
	int w = (img.Width / init_sample);
	int h = (img.Height / init_sample);
	int s = (init_sample);

	// Calculate approximated determinant of hessian values
	if (octaves >= 1)
	{
		responseMap->push_back(ResponseLayer(w, h, s,  9));
		responseMap->push_back(ResponseLayer(w, h, s, 15));
		responseMap->push_back(ResponseLayer(w, h, s, 21));
		responseMap->push_back(ResponseLayer(w, h, s, 27));
	}
	 
	if (octaves >= 2)
	{
		responseMap->push_back(ResponseLayer(w / 2, h / 2, s * 2, 39));
		responseMap->push_back(ResponseLayer(w / 2, h / 2, s * 2, 51));
	}

	if (octaves >= 3)
	{
		responseMap->push_back(ResponseLayer(w / 4, h / 4, s * 4, 75));
		responseMap->push_back(ResponseLayer(w / 4, h / 4, s * 4, 99));
	}

	if (octaves >= 4)
	{
		responseMap->push_back(ResponseLayer(w / 8, h / 8, s * 8, 147));
		responseMap->push_back(ResponseLayer(w / 8, h / 8, s * 8, 195));
	}

	if (octaves >= 5)
	{
		responseMap->push_back(ResponseLayer(w / 16, h / 16, s * 16, 291));
		responseMap->push_back(ResponseLayer(w / 16, h / 16, s * 16, 387));
	}

	  // Extract responses from the image
	for (unsigned int i = 0; i < responseMap->size(); ++i)
	{
		buildResponseLayer((*responseMap)[i]);
	}
}


/// <summary>
/// Build Responses for a given ResponseLayer
/// </summary>
/// <param name="rl"></param>
void FastHessian::buildResponseLayer(ResponseLayer rl)
{
	int step = rl.step;                      // step size for this filter
	int b = (rl.filter - 1) / 2;             // border for this filter
	int l = rl.filter / 3;                   // lobe for this filter (filter size / 3)
	int w = rl.filter;                       // filter size
	float inverse_area = (float)1.0 / (w * w);       // normalisation factor
	float Dxx, Dyy, Dxy;

	for (int r, c, ar = 0, index = 0; ar < rl.height; ++ar)
	{
		for (int ac = 0; ac < rl.width; ++ac, index++)
		{
			// get the image coordinates
			r = ar * step;
			c = ac * step;

			// Compute response components
			Dxx = img.BoxIntegral(r - l + 1, c - b, 2 * l - 1, w)
			    - img.BoxIntegral(r - l + 1, c - l / 2, 2 * l - 1, l) * 3;
			Dyy = img.BoxIntegral(r - b, c - l + 1, w, 2 * l - 1)
			    - img.BoxIntegral(r - l / 2, c - l + 1, l, 2 * l - 1) * 3;
			Dxy = + img.BoxIntegral(r - l, c + 1, l, l)
			      + img.BoxIntegral(r + 1, c - l, l, l)
			      - img.BoxIntegral(r - l, c - l, l, l)
			      - img.BoxIntegral(r + 1, c + 1, l, l);

			// Normalise the filter responses with respect to their size
			Dxx *= inverse_area;
			Dyy *= inverse_area;
			Dxy *= inverse_area;

			// Get the determinant of hessian response & laplacian sign
			rl.responses[index] = (Dxx * Dyy - (float)0.81 * Dxy * Dxy);
			rl.laplacian[index] = (char)(Dxx + Dyy >= 0 ? 1 : 0);
		}
	}
}


/// <summary>
/// Test whether the point r,c in the middle layer is extremum in 3x3x3 neighbourhood
/// </summary>
/// <param name="r">Row to be tested</param>
/// <param name="c">Column to be tested</param>
/// <param name="t">Top ReponseLayer</param>
/// <param name="m">Middle ReponseLayer</param>
/// <param name="b">Bottome ReponseLayer</param>
/// <returns></returns>
bool FastHessian::isExtremum(int r, int c, ResponseLayer t, ResponseLayer m, ResponseLayer b)
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


/// <summary>
/// Interpolate scale-space extrema to subpixel accuracy to form an image feature
/// </summary>
/// <param name="r"></param>
/// <param name="c"></param>
/// <param name="t"></param>
/// <param name="m"></param>
/// <param name="b"></param>
void FastHessian::interpolateExtremum(int r, int c, ResponseLayer t, ResponseLayer m, ResponseLayer b)
{
	double * D    = BuildDerivative(r, c, t, m, b);
	double ** H  = BuildHessian(r, c, t, m, b);
	double ** Hi = Inverse(H);

	if (Hi != NULL) {
		double * Of = MMM_neg_3x3_3x1(Hi, D);

		// get the offsets from the interpolation
		double O[3] = { Of[0], Of[1], Of[2] };

		// get the step distance between filters
		int filterStep = (m.filter - b.filter);
	 
		// If point is sufficiently close to the actual extremum
		if (abs(O[0]) < (float)0.5 && abs(O[1]) < (float)0.5 && abs(O[2]) < (float)0.5)
		{
			IPoint ipt;
			ipt.x = (float)((c + O[0]) * t.step);
			ipt.y = (float)((r + O[1]) * t.step);
			ipt.scale = (float)((0.1333f) * (m.filter + O[2] * filterStep));
			ipt.laplacian = (int)(m.getLaplacian(r,c,t));
			ipts->push_back(ipt);
		}

		delete[] Of;
	}


	delete[] H[0];
	delete[] H[1];
	delete[] H[2];
	delete[] Hi[0];
	delete[] Hi[1];
	delete[] Hi[2];
	delete[] H;
	delete[] Hi;

}

/// <summary>
/// Build Matrix of First Order Scale-Space derivatives
/// </summary>
/// <param name="octave"></param>
/// <param name="interval"></param>
/// <param name="row"></param>
/// <param name="column"></param>
/// <returns>3x1 Matrix of Derivatives</returns>
double * FastHessian::BuildDerivative(int r, int c, ResponseLayer t, ResponseLayer m, ResponseLayer b)
{
	double dx, dy, ds;

	dx = (m.getResponse(r, c + 1, t) - m.getResponse(r, c - 1, t)) / (float)(2.0);
	dy = (m.getResponse(r + 1, c, t) - m.getResponse(r - 1, c, t)) / (float)(2.0);
	ds = (t.getResponse(r, c) - b.getResponse(r, c, t)) / (float)(2.0);

	double * D = new double[3];
	D[0] = dx;
	D[1] = dy;
	D[2] = ds;
	return D;
}


/// <summary>
/// Build Hessian Matrix 
/// </summary>
/// <param name="octave"></param>
/// <param name="interval"></param>
/// <param name="row"></param>
/// <param name="column"></param>
/// <returns>3x3 Matrix of Second Order Derivatives</returns>
double ** FastHessian::BuildHessian(int r, int c, ResponseLayer t, ResponseLayer m, ResponseLayer b)
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

	double ** H = new double * [3];
	H[0] = new double[3];
	H[1] = new double[3];
	H[2] = new double[3];

	H[0][0] = dxx;
	H[0][1] = dxy;
	H[0][2] = dxs;
	H[1][0] = dxy;
	H[1][1] = dyy;
	H[1][2] = dys;
	H[2][0] = dxs;
	H[2][1] = dys;
	H[2][2] = dss;
	return H;
}

double ** FastHessian::Inverse(double ** m)
{
	double det = ((m[0][0]*m[1][1]*m[2][2]) + (m[0][1]*m[1][2]*m[2][0]) + (m[0][2]*m[1][0]*m[2][1]) -
				  (m[0][2]*m[1][1]*m[2][0]) - (m[0][1]*m[1][0]*m[2][2]) - (m[0][0]*m[1][2]*m[2][1]));

	if (det == 0) return NULL;

	double ** invert = new double * [3];
	invert[0] = new double[3];
	invert[1] = new double[3];
	invert[2] = new double[3];
	double A =      (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
	double B = -1 * (m[1][0] * m[2][2] - m[2][0] * m[1][2]) / det;
	double C =      (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det;
	double D = -1 * (m[0][1] * m[2][2] - m[0][2] * m[2][1]) / det;
	double E =      (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
	double F = -1 * (m[0][0] * m[2][1] - m[0][1] * m[2][0]) / det;
	double G =      (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
	double H = -1 * (m[0][0] * m[1][2] - m[0][2] * m[1][0]) / det;
	double K =      (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
	
	invert[0][0] = A;
	invert[0][1] = D;
	invert[0][2] = G;
	invert[1][0] = B;
	invert[1][1] = E;
	invert[1][2] = H;
	invert[2][0] = C;
	invert[2][1] = F;
	invert[2][2] = K;

	return invert;

}

double * FastHessian::MMM_neg_3x3_3x1(double ** A, double * B)
{
	double * C = new double[3];
	double a = -1 * (A[0][0] * B[0] + A[0][1] * B[1] + A[0][2] * B[2]);
	double b = -1 * (A[1][0] * B[0] + A[1][1] * B[1] + A[1][2] * B[2]);
	double c = -1 * (A[2][0] * B[0] + A[2][1] * B[1] + A[2][2] * B[2]);

	C[0] = a;
	C[1] = b;
	C[2] = c;
	
	return C;
}
