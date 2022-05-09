
#include <complex>
#include <cmath>

	// Sinlike
	// a = power of since argument - acts like the number of regions where points accumulate periodically
	// n = nodes
	// i = [1,n]

	static double Sinlike(int i, int n, double a, double rmax, double rmin) {
		return rmin+rmax*sin(M_PI/2. * pow(double(i-1)/double(n-1), a));
	}
	// Exponentential
	// becomes linear as G->0. Points accumulate at rmin as g->infty
	// n = nodes
	// i = [1,n]
	// g = accumulation factor
	static double Exponential(int i, int n, double g, double rmin, double rmax) {
		return rmin + (rmax - rmin)*(exp(g*(i-1)/(n-1)) - 1.)/(exp(g) - 1.);
	}


	// Linear Parabolic
	// quadratic near rmin, linear at large distances
	// i0 = break point
	// n = nodes
	// i = [1,n]
	static double LinearParabolic(int i, int n, int i0, double rmin, double rmax) {
		double r0 = (rmax*(i0 - 1) + rmin*(n-i0)) / (2*n-i0 - 1);			// the node at i0
		double a = (r0 - rmin)/(i0-1)/(i0-1);
		double b = (rmax - r0)/(n-i0);


		if (i < i0) 
			return rmin + a*(i-1)*(i-1);
		return r0 + b*(i-i0);
	}
	// linear
	static double Linear(int i, int n, double rmin, double rmax) {
		return rmin + (rmax - rmin)/double(n-1)*(i-1);
	}
