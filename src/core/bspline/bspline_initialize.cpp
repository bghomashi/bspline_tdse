#include "bspline.h"
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std::complex_literals;

namespace Basis {

	int BSpline::Initialize(int order, int nodes, double xmin, double xmax, Sequence seq, const ECS& ecs, double param) {
		_nodes = nodes;
		_order = order;
		_numBSplines = (nodes-1) + (order-1);
		_skipFirst = false;
		_skipLast = false;
		_grid = std::vector<double>(nodes);
		_ecs_grid = std::vector<complex>(nodes);


		// initialize grid -- UGLY. MOVE THIS OUT OF THIS FUNCTION
		if (seq == Linear) {
			for (int i = 0; i < nodes; i++)
				_grid[i] = xmin + ((xmax - xmin) * (double)i) / double(nodes - 1.);
		} else if (seq == Exponential) {
			double g = param;											// parameter
			for (int i = 0; i < nodes; i++)
				_grid[i] = xmin + (xmax - xmin)*(exp(g*i/(nodes-1)) - 1.)/(exp(g) - 1.);
		} else if (seq == ParabolicLinear) {
			double x0 = param;
			double i0 = std::floor(2.*(nodes-1) / (1.+(xmax - xmin)/(x0 - xmin)));
			// double i0 = std::min(30, nodes);						// parameter
			double a0 = xmax / i0 /(2.*(nodes-1) - i0);
			double a1 = 2.*xmax /(2.*(nodes-1) - i0);
			double a2 = -xmax*i0 /(2.*(nodes-1) - i0);

			for (int i = 0; i < nodes; i++) {
				if (i < i0) 
					_grid[i] = xmin + a0*i*i;
				else
					_grid[i] = xmin + a2 + a1*i;
			}
		} else if (seq == Sinlike) {
			// not tested
			double a = param;											// parameter
			for (int i = 0; i < nodes; i++)
				_grid[i] = xmin + xmax*sin(M_PI/2. * pow(double(i)/double(nodes-1), a));
		// std::ofstream file("data/grid_sinlike.dat");
		// for (int i = 0; i < nodes; i++)
		// 	file << i << "\t" << _grid[i] << "\n";
		// file.close();
		// exit(0);
		}


        // initialize ecs
		assert((ecs.r0 > 0. && ecs.r0 <= 1.) && "ecs must be between 0 and 1.");
		if (ecs.r0 == 1.) {
			_ecs.r0 = xmax;
			_ecs.theta = 0.0;				// no ecs
		} else {
			_ecs.theta = ecs.theta;
			_ecs.r0 = _grid[whichInterval(ecs.r0*(xmax-xmin)+xmin)];	// clamp ecs to a knot
		}


		for (int i = 0; i < _grid.size(); i++)
			_ecs_grid[i] = _ecs.R(_grid[i]);



		// initialize knots
		_knots = std::vector<complex>(nodes + 2*(order-1));
		for (int i = 0; i < order; i++)
			_knots[i] = _ecs.R(_grid[0]);
		for (int i = order; i < nodes + order - 2; i++)
			_knots[i] = _ecs.R(_grid[i - order + 1]);
		for (int i = nodes + order - 2; i < nodes + 2*(order-1); i++)
			_knots[i] = _ecs.R(_grid[nodes - 1]);



		// initialize Gauss-Legendre quadrature
		_glQuad.Initialize();

		InitializeFactorials();
		InitializeBSplines();

		return 0;
	}

	void BSpline::InitializeFactorials() {
		// initialize factorials
		_factorial.resize(DIMFACT + 1);
		_partialFactorial.resize((DIMFACT + 1) * (DIMFACT + 1));

		_factorial[0] = 1.;
		for (int i = 1; i <= DIMFACT; i++)
			_factorial[i] = _factorial[i - 1] * (double)i;
		// this wastes about 50% of the memory: 0.5*127^2*sizeof(double) ~ 64 kb
		// initialize partial factorial array
		for (int i = 0; i <= DIMFACT; i++)
			for (int j = 0; j <= DIMFACT; j++)
				_partialFactorial[i + j * DIMFACT] = 1.;


		for (int n = 1; n <= DIMFACT; n++)
			for (int k = 1; k <= n; k++)
				for (int i = n - k + 1; i <= n; i++)
					_partialFactorial[n + k * DIMFACT] *= (double)i;
	}
	void BSpline::InitializeBSplines() {
		_bsCoeffs.resize(_numBSplines, BSplineCoeff(_order-1));

		std::vector<complex> vec(_order + 1);

		for (int bs = 0; bs < _numBSplines; bs++) {
			for (int i = 0; i <= _order; i++)
				vec[i] = _knots[bs + i];

			InitializeBSpline(vec, _bsCoeffs[bs]);
		}
	}
	void BSpline::InitializeBSpline(const std::vector<complex>& knots, BSplineCoeff& out) {
		CoeffMatrix coeff(_order-1, knots.size());

		// For 1st order, there is only one interval in the bspline and it is constant.
		for (int i = 0; i < _order; i++) {
			auto& bs = coeff(0, i);
			bs(0, 0) = 1.;
		}

		// For each order after 1
		for (int k = 2; k <= _order /*-1??*/; k++) {
			InitializeBSplineOfOrder(knots, k, coeff);
		}
		auto& bs = coeff(_order-1, 0);
		for (int interval = 0; interval < _order; interval++)
			for (int d = 0; d < _order; d++)
				out(interval, d) = bs(interval, d);
	}

	void BSpline::InitializeBSplineOfOrder(const std::vector<complex>& knots, int order, CoeffMatrix& coeff) {
		// For each BSpline of this order (num_knots - order - 1)
		for (int tbs = 0; tbs < _order - order + 1; tbs++) {
			PropagateCoeffients(knots, tbs, order, coeff);
		}
	}

	void BSpline::PropagateCoeffients(const std::vector<complex>& knots, int tbs, int order, CoeffMatrix& coeff) {
		auto& bs = coeff(order-1, tbs);
		complex tbsSpan, A;
		// Each BSpline has is a combination of two lower order bsplines
		// first the left most
		tbsSpan = knots[tbs + order-1] - knots[tbs];	//total span of sub-BSpline
		if (std::abs(tbsSpan) > NOD_THRESHOLD) {
			const auto& bs1 = coeff(order - 2, tbs);		// get the list of coefficients for the left most BSpline
			// now we need to go through each interval of that sub-BSpline
			for (int interval = 0; interval < order - 1; interval++) {

				// ---- constant component ----
				A = (knots[tbs + interval] - knots[tbs]) / tbsSpan;		// normalize this interval
				for (int d = 0; d < order-1; d++)
					bs(interval, d) += bs1(interval, d) * A;
				// ---- monomial component ----
				A = (knots[tbs + interval + 1] - knots[tbs + interval]) / tbsSpan;		// normalize this interval
				for (int d = 0; d < order-1; d++)
					bs(interval, d + 1) += bs1(interval, d) * A;
			}

		}
		// now the right bspline
		tbsSpan = knots[tbs + order] - knots[tbs + 1];
		if (std::abs(tbsSpan) > NOD_THRESHOLD) {
			const auto& bs2 = coeff(order - 2, tbs + 1);		// get the list of coefficients for the left most BSpline
			// now we need to go through each interval of that sub-BSpline
			for (int interval = 0; interval < order - 1; interval++) {

				// ---- constant component ----
				A = (knots[tbs + order] - knots[tbs + interval + 1]) / tbsSpan;		// normalize this interval
				for (int d = 0; d < order-1; d++)
					bs(interval + 1, d) += bs2(interval, d) * A;
				// ---- monomial component ----
				A = -(knots[tbs + interval + 2] - knots[tbs + interval + 1]) / tbsSpan;		// normalize this interval
				for (int d = 0; d < order-1; d++)
					bs(interval + 1, d + 1) += bs2(interval, d) * A;
			}
		}
	}
}