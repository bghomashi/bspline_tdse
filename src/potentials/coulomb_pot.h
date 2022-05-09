#pragma once

#include "tdse/potential.h"

class CoulombPotential : public Potential {
    double _Z;
    double _x, _y, _z;
public:
    CoulombPotential();
    CoulombPotential(double Z);
    CoulombPotential(double Z, double x, double y, double z);

    // FIX: maybe this way of thinking about it is not very efficient
    void FillMatrix(const Basis::BSpline& basis, Matrix m, int N, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradX(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradY(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradZ(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
};