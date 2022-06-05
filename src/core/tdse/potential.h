#pragma once

#include "maths/maths.h"
#include "bspline/bspline.h"
#include <cassert>
#include <memory>

class Potential {
protected:
    bool _isCentral;
    bool _isAxial;
public:
    typedef std::shared_ptr<Potential> Ptr_t;

    Potential() : _isCentral(true), _isAxial(true) {};

    
    inline bool isCentral() const {
        return _isCentral;
    }
    inline bool isAxial() const {
        return _isAxial;
    }

    virtual double operator() (double x, double y, double z) const = 0;
    virtual void FillMatrix(const Basis::BSpline& basis, Matrix m, int N, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0}) = 0;
    virtual void FillMatrixGradX(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0}) = 0;
    virtual void FillMatrixGradY(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0}) = 0;
    virtual void FillMatrixGradZ(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0}) = 0;
};