#include "bspline.h"

namespace Basis {
	BSpline::BSplineCoeff::BSplineCoeff(int max_order) {
		d.resize(max_order + 1);			// intervals
		for (auto& c : d)
			c.resize(max_order + 1, 0);		//degree coeffs
	}

	complex BSpline::BSplineCoeff::operator() (int interval, int degree) const {
		if (interval >= d.size()) return 0.;
		return d[interval][degree];
	}
	complex& BSpline::BSplineCoeff::operator() (int interval, int degree) {
		return d[interval][degree];
	}
	BSpline::CoeffMatrix::CoeffMatrix(size_t max_order, size_t max_splines) {
		d.resize(max_order + 1);									// each order [0 ... max]
		for (int o = 0; o < max_order + 1; o++) {
			d[o].resize(max_splines - o - 1, BSplineCoeff(o));		// has this many
		}
	}

	BSpline::BSplineCoeff& BSpline::CoeffMatrix::operator() (size_t order, size_t spline) {
		return d[order][spline];
	}
}