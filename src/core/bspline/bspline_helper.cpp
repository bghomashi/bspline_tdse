#include "bspline.h"


namespace Basis {
	int BSpline::getNumBSplines() const {
		int num = _numBSplines;
		if (_skipFirst) num--;
		if (_skipLast) num--;
		return num;
	}
	int BSpline::getOrder() const {
		return _order;
	}
	int BSpline::whichInterval(double x) const {
		for (int i = 0; i < (int)_grid.size(); i++) {
			if (_grid[i] <= x && x <= _grid[i + 1])
				return i;
		}
		return -1;
	}
	
	const std::vector<double>& BSpline::getGrid() const {
		return _grid;
	}
	void BSpline::setSkipFirst(bool flag) {
		_skipFirst = flag;
	}
	void BSpline::setSkipLast(bool flag) {
		_skipLast = flag;
	}
}