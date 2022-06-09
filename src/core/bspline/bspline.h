#pragma once

#include "maths/maths.h"
#include <functional>
#include <vector>

namespace Basis {  
    static constexpr double NOD_THRESHOLD = 1.E-15;								// grid space that is considered zero
    static constexpr complex I =  complex(0,1);

	class GaussQuadrature {
		static constexpr double GAUSS_POINTS_CONVERGENCE_THRESHOLD = 2.E-15;	
		const int NGAUSS = 64;													// how many points

		std::vector<double> _points, _weights;
	public:
		bool Initialize();
		complex Integrate(complex xmin, complex xmax, std::function<complex(complex)> f) const;
	};

    class BSpline {
    public:
		enum Sequence {
			Linear,
			Exponential,
			Sinlike,
			ParabolicLinear
		};
		struct ECS {
			double r0;					// [0, 1]
			double theta;				// radians 

			complex q(double x) const {
				if (x <= r0)
					return 1;
				return std::exp(I*theta);
			}
			complex R(double x) const {
				if (x <= r0)
					return x;
				return r0 + (x-r0)*std::exp(I*theta);
			}
			double x(complex R) const {
				if (std::imag(R) == 0)
					return std::real(R);
				return r0 + std::real((R - r0)*std::exp(-I*theta));
			}
		};
    private:
		struct BSplineCoeff {
			std::vector<std::vector<complex>> d;

			BSplineCoeff(int max_order);
			complex operator() (int interval, int degree) const;
			complex& operator() (int interval, int degree);
		};
		struct CoeffMatrix {
			std::vector<std::vector<BSplineCoeff>> d;
			CoeffMatrix(size_t order, size_t num_knots);
			BSplineCoeff& operator() (size_t order, size_t bspline);
		};

        GaussQuadrature _glQuad;
        ECS _ecs;

        
		std::vector<double> _grid;                  // real locations of the nodes
		std::vector<complex> _ecs_grid, _knots;     // complex nodes, knots
		
		int _nodes;                                 // number of nodes
		int _order;                                 // order of bsplines
		int _numBSplines;                           // number of bsplines
		bool _skipFirst, _skipLast;                 // assert 0 boundary by excluding first/last bspline

        
		std::vector<BSplineCoeff> _bsCoeffs;        // polynomial coefficients for each interval
        
        // factorial storage
		const int DIMFACT = 127;
		std::vector<double> _partialFactorial, _factorial;
        
        
		void InitializeFactorials();
		void InitializeBSplines();
		void InitializeBSpline(const std::vector<complex>& knots, BSplineCoeff& out);
		void InitializeBSplineOfOrder(const std::vector<complex>& knots, int order, CoeffMatrix& coeff);
		void PropagateCoeffients(const std::vector<complex>& knots, int tbs, int order, CoeffMatrix& coeff);
    public:
		int Initialize(int order, int nodes, double xmin, double xmax, Sequence seq = Linear, const ECS& ecs = {0.9, Pi/4.}, double param = 0.);

		// get the bs'th Bspline (derivative dn)
		complex bspline(double x, int bs, int dn = 0) const;                // x-before ecs rotation
		complex bspline(complex x, int bs, int dn = 0) const;               // x-after ecs rotation
		std::vector<complex> getBSpline(const std::vector<double>& x, int bs, int dn = 0) const;// x-before ecs rotation

		// evaluate function expanded in Bspline basis at x - dn is the order of the derivative
        // x-before ecs rotation
		complex FunctionEvaluate(double x, const std::vector<complex>& fc, int dn = 0) const;
		std::vector<complex> FunctionEvaluate(const std::vector<double>& x, const std::vector<complex>& fc, int dn = 0) const;

		// integrate a function over two Bsplines over the whole grid
        // the parameter passed to 'f' is (complex) x
		complex Integrate(int bs1, int bs2, int dn1 = 0, int dn2 = 0) const;
		complex Integrate(int bs1, int bs2, std::function<complex(complex)> f, int dn1 = 0, int dn2 = 0) const;
		
        // integrate a function over restricted bounds
		complex Integrate(double xmin, double xmax, int bs1, int bs2, std::function<complex(complex)> f, int dn1 = 0, int dn2 = 0) const;

		// helper functions
		const std::vector<double>& getGrid() const;
		int getNumBSplines() const;
		int getOrder() const;
		int whichInterval(double x) const;
		void setSkipFirst(bool flag = true);
		void setSkipLast(bool flag = true);
    };   
}