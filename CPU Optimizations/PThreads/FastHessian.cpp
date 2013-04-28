#include "IPoint.h"
#include <cmath>
#include <vector>
#include "FastHessian.h"
#include "IntegralImage.h"
#include "parameters.h"
#include <pthread.h>

struct thread_data {
	int thread_id;
	FastHessian::ResponseLayer * rl;
	IntegralImage * img;
	int step;
	int b;
	int l;
	int w;
	float inverse_area;
};


/// <summary>
/// Static one-call do it all method
/// </summary>
/// <param name="thresh"></param>
/// <param name="octaves"></param>
/// <param name="init_sample"></param>
/// <param name="img"></param>
/// <returns></returns>
std::vector<IPoint> * FastHessian::getIpoints(float thresh, int octaves, int init_sample, IntegralImage * img)
{
	FastHessian fh(thresh, octaves, init_sample, img);
	return fh.getIpoints();
}


/// <summary>
/// Constructor with parameters
/// </summary>
/// <param name="thresh"></param>
/// <param name="octaves"></param>
/// <param name="init_sample"></param>
/// <param name="img"></param>
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


/// <summary>
/// Build Responses for a given ResponseLayer
/// </summary>
/// <param name="rl"></param>
void FastHessian::buildResponseLayer(ResponseLayer &rl)
{
	void * buildResponseLayer_work(void * threadarg);
	int step = rl.step;                      // step size for this filter
	int b = (rl.filter - 1) / 2;             // border for this filter
	int l = rl.filter / 3;                   // lobe for this filter (filter size / 3)
	int w = rl.filter;                       // filter size
	float inverse_area = (float)1.0 / (w * w);       // normalisation factor

	pthread_t threads[NUM_THREADS];
	struct thread_data data[NUM_THREADS];

	for (int t = 0; t < NUM_THREADS; t++)
	{
		data[t].thread_id = t;
		data[t].rl = &rl;
		data[t].img = img;
		data[t].step = step;
		data[t].b = b;
		data[t].l = l;
		data[t].w = w;
		data[t].inverse_area = inverse_area;

		int rc = pthread_create(&threads[t], NULL, buildResponseLayer_work, (void*) &(data[t]));
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
}

void * buildResponseLayer_work(void * threadarg)
{
	struct thread_data * my_data = (struct thread_data *) threadarg;
	int thread_id = my_data->thread_id;
	FastHessian::ResponseLayer * rl = my_data->rl;
	IntegralImage * img = my_data->img;
	int step = my_data->step;                      // step size for this filter
	int b = my_data->b;             // border for this filter
	int l = my_data->l;                   // lobe for this filter (filter size / 3)
	int w = my_data->w;                       // filter size
	float inverse_area = my_data->inverse_area;       // normalisation factor
	float Dxx, Dyy, Dxy;

	unsigned int start = (rl->height / NUM_THREADS) * thread_id;
	unsigned int end = (rl->height / NUM_THREADS) + start;
	if (end > rl->height) end  = rl->height;

	for (int r, c, ar = start, index = start * rl->width; ar < end; ++ar)
	{
		for (int ac = 0; ac < rl->width; ++ac, index++)
		{

			// get the image coordinates
			r = ar * step;
			c = ac * step;

			// Compute response components
			Dxx = img->BoxIntegral(r - l + 1, c - b, 2 * l - 1, w)
			    - img->BoxIntegral(r - l + 1, c - l / 2, 2 * l - 1, l) * 3;
			Dyy = img->BoxIntegral(r - b, c - l + 1, w, 2 * l - 1)
			    - img->BoxIntegral(r - l / 2, c - l + 1, l, 2 * l - 1) * 3;
			Dxy = + img->BoxIntegral(r - l, c + 1, l, l)
			      + img->BoxIntegral(r + 1, c - l, l, l)
			      - img->BoxIntegral(r - l, c - l, l, l)
			      - img->BoxIntegral(r + 1, c + 1, l, l);

			// Normalise the filter responses with respect to their size
			Dxx *= inverse_area;
			Dyy *= inverse_area;
			Dxy *= inverse_area;

			// Get the determinant of hessian response & laplacian sign
			rl->responses[index] = (Dxx * Dyy - (float)0.81 * Dxy * Dxy);
			rl->laplacian[index] = (unsigned char)(Dxx + Dyy >= 0 ? 1 : 0);
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


/// <summary>
/// Interpolate scale-space extrema to subpixel accuracy to form an image feature
/// </summary>
/// <param name="r"></param>
/// <param name="c"></param>
/// <param name="t"></param>
/// <param name="m"></param>
/// <param name="b"></param>
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

/// <summary>
/// Build Matrix of First Order Scale-Space derivatives
/// </summary>
/// <param name="octave"></param>
/// <param name="interval"></param>
/// <param name="row"></param>
/// <param name="column"></param>
/// <returns>3x1 Matrix of Derivatives</returns>
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


/// <summary>
/// Build Hessian Matrix 
/// </summary>
/// <param name="octave"></param>
/// <param name="interval"></param>
/// <param name="row"></param>
/// <param name="column"></param>
/// <returns>3x3 Matrix of Second Order Derivatives</returns>
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
