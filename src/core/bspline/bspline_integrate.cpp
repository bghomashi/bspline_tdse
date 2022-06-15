#include "bspline.h"
#include <iostream>
#include <cmath>

namespace Basis {
	complex BSpline::Integrate(int bs1, int bs2, int dn1, int dn2) const {
		return Integrate(bs1, bs2, [](complex r) -> complex {
				return 1.;
			}, dn1, dn2);
	}
	complex BSpline::Integrate(int bs1, int bs2, std::function<complex(complex)> f, int dn1, int dn2) const {
		double xmin = _grid.front(), xmax = _grid.back();
		return Integrate(xmin, xmax, bs1, bs2, f, dn1, dn2);
	}


/* --- this is probably pointless? 
TODO: check

function KahanSum(input)
    var sum = 0.0                    // Prepare the accumulator.
    var c = 0.0                      // A running compensation for lost low-order bits.

    for i = 1 to input.length do     // The array input has elements indexed input[1] to input[input.length].
        var y = input[i] - c         // c is zero the first time around.
        var t = sum + y              // Alas, sum is big, y small, so low-order digits of y are lost.
        c = (t - sum) - y            // (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
        sum = t                      // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
    next i                           // Next time around, the lost low part will be added to y in a fresh attempt.

    return sum
*/
	complex BSpline::Integrate(double xmin, double xmax, int bs1, int bs2, std::function<complex(complex)> f, int dn1, int dn2) const {
		// easy out - no overlap of the bsplines
		// std::cout << "bs1 " << bs1 << " bs2 " << bs2 << std::endl;
		// std::cout << "bs1 - _order + 1 " << bs1 - _order + 1 << " bs2 - _order + 1 " << bs2 - _order + 1 << std::endl;

		if (bs1 - _order + 1 > bs2 || bs2 - _order + 1 > bs1)return 0;
		
		// first find the intervals we need to integrate over. This is the same on any reasonable ecs contour.
		int IntervalMin = whichInterval(xmin);
		int IntervalMax = whichInterval(xmax);
		IntervalMin = std::max<int>(IntervalMin, std::max<int>(bs1 - _order + 1, bs2 - _order + 1));
		IntervalMax = std::min<int>(IntervalMax, std::min<int>(bs1, bs2));

		// we are actuall 

		// here i used some fancy sum to account for round off error... not sure if it matters		
		complex total = 0., elem = 0., y = 0., t = 0., c = 0.;
		for (int i = IntervalMin; i <= IntervalMax; i++) {	// each interval
			// evaluate the integral on the (complex) contour
			complex lower_bound = _ecs.R(std::max(xmin, _grid[i]));		// start from xmin
			complex upper_bound = _ecs.R(std::min(xmax, _grid[i + 1]));	// stop at xmax

			elem = _glQuad.Integrate(lower_bound, upper_bound, [=](complex x) {
				return f(x) * (bspline(x, bs1, dn1) * bspline(x, bs2, dn2));
			});
			
			y = elem - c;
			t = total + y;
			c = (t - total) - y;
			total = t;
		}
		return total;
	}
}