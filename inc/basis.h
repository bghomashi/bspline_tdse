#pragma once

#include <vector>
#include <cassert>

#include <functional>

#include "matrix.h"
#include "maths.h"

class Basis {

public:

};

class FiniteDifference : public Basis {
    int _dimensions;
    double _step;
    std::vector<int> _numPoints;
    std::vector<double> _xmin, _xmax;
public:
    FiniteDifference(double xmin, double xmax, double step) : 
        _step(step), _dimensions(1), _numPoints(_dimensions), _xmin(_dimensions), _xmax(_dimensions) {
        _xmin[0] = xmin;
        _xmax[0] = xmax;
        _numPoints[0] = int((xmax - xmin) / step + 1);
    }
    FiniteDifference(const std::vector<double>& xmin, const std::vector<double>& xmax, const double step) :
        _step(step), _dimensions(xmin.size()), _numPoints(_dimensions), _xmin(_dimensions), _xmax(_dimensions) {
        assert(xmin.size() == xmax.size() && "xmin and xmax must be the same length.");
        for (int id : range(0, _dimensions)) {
            _xmin[id] = xmin[id];
            _xmax[id] = xmax[id];
            _numPoints[id] = int((xmax[id] - xmin[id])/step + 1);
        }
    }

    int GetNumNodes(int dim) const {
        return _numPoints[dim];
    }

    void FillLaplacian(Matrix& out /*, int order = 2*/) {
        // total size of matrix numPoints^dim
        int totalNumPoints = 1;
        for (int id : range(0, _dimensions))
            totalNumPoints *= _numPoints[id];

        std::vector<int> ix(_dimensions);
        for (int ip : range(0, totalNumPoints)) {
            // ip = x + Nx*(y + Ny*z)

            int temp_ip = ip;
            for (int id : range(0, _dimensions)) {
                ix[id] = temp_ip % _numPoints[id];
                temp_ip /= _numPoints[id];
            }

            int dimScale = 1;                    // this will keep track of which matrix block we are looking at
            for (int id : range(0, _dimensions)) {
                if (ix[id]-1 >= 0)               // look back one point along this dimension
                    out.Set(ip, ip-dimScale, 1.);
                out.Add(ip, ip, -2.);                   // always add to diagonal
                if (ix[id]+1 < _numPoints[id])   // look forward one point along this dimension
                    out.Set(ip, ip+dimScale, 1.);

                dimScale *= _numPoints[id];      // now skip over this whole dimension for the next points
            }
        }
    }
    void FillFunction(std::function<complex(const std::vector<complex>&)> f, Matrix& out) {
        // total size of matrix numPoints^dim
        int totalNumPoints = 1;
        for (int id : range(0, _dimensions))
            totalNumPoints *= _numPoints[id];

        std::vector<int> ix(_dimensions);
        for (int ip : range(0, totalNumPoints)) {
            // ip = x + Nx*(y + Ny*z)

            int temp_ip = ip;
            for (int id : range(0, _dimensions)) {
                ix[id] = temp_ip % _numPoints[id];
                temp_ip /= _numPoints[id];
            }

            std::vector<complex> r(ix.size());
            for (int id : range(0, _dimensions))
                r[id] = ix[id]*_step + _xmin[id];
                
            out.Add(ip, ip, f(r));                   // always add to diagonal
        }
    }
};

class BSpline : public Basis{
public:
};

