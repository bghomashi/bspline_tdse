#pragma once

#include "tdse/potential.h"

class ExponentialPotential : public Potential {
    double _Z, _D;
    double _x, _y, _z;
public:
    ExponentialPotential();
    ExponentialPotential(double Z, double decay);
    ExponentialPotential(double Z, double decay, double x, double y, double z);

    double operator() (double x, double y, double z) const;
    void FillMatrix(const Basis::BSpline& basis, Matrix m, int N, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradX(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradY(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradZ(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
};