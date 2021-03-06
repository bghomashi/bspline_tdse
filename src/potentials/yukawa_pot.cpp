#include "potentials/yukawa_pot.h"
#include "utility/index_manip.h"
#include "utility/logger.h"
#include "utility/spherical_harmonics.h"

#include <iostream>
#include <functional>
#include <string>


static void BuildExpInvR(const Basis::BSpline& basis, int N, double Z, double D, std::vector<complex>& invR);
static void BuildExpInvRR(const Basis::BSpline& basis, int N, double Z, double D, std::vector<complex>& invRR);
static void FillBlock( int l1, int m1, int l2, int m2,
                int lmax, int mmax, int N, int order,
                const std::vector<int>& Ms,
                const std::vector<int>& mRows,
                const std::vector<complex>& invR,
                const std::vector<complex>& invRR,
                std::function<complex(int, int, int, int)> YlmYlm,
                Matrix m);



YukawaPotential::YukawaPotential() : _Z(0), _D(0), _x(0), _y(0), _z(0) {}
YukawaPotential::YukawaPotential(double Z, double decay) : _Z(Z), _D(decay), _x(0), _y(0), _z(0) {}
YukawaPotential::YukawaPotential(double Z, double decay, double x, double y, double z) :
    _Z(Z), _D(decay), _x(x), _y(y), _z(z) {
    if (z != 0) {
        _isCentral = false;
        if (x != 0 || y != 0)
            _isAxial = false;
    }
}
double YukawaPotential::operator() (double x, double y, double z) const {
    double r = std::sqrt (x*x + y+y + z*z);
    return (-_Z/r)*std::exp(-_D*r);
}

void YukawaPotential::FillMatrix(const Basis::BSpline& basis, Matrix m, int N, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = basis.getOrder();

    if (_isCentral) {
        // the same for each L-M block so cache it
        std::vector<complex> invR(N*N);
        BuildExpInvR(basis, N, _Z, _D, invR);

        m->FillBandedBlock(order-1, N, [&](int row, int col) {
            int i, j, l1, l2, m1, m2;
            ILMFrom(row, i, l1, m1, N, Ms, mRows);
            ILMFrom(col, j, l2, m2, N, Ms, mRows);
            return invR[i + j*N];
        });
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
}
void YukawaPotential::FillMatrixGradX(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = basis.getOrder();
    int mmax = Ms.back();

    if (_isCentral) {
        // the same for each L-M block so cache it
        std::vector<complex> invR(N*N), invRR(N*N);
        BuildExpInvR(basis, N, _Z, _D, invR);
        for (auto& i : invR) i*=-_D;                 // from product rule
        BuildExpInvRR(basis, N, _Z, _D, invRR);

        // for each m-block (block rows)
        for (int m1 : Ms) {
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
                FillBlock(l1, m1, l1+1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmXYlm, m);
                FillBlock(l1, m1, l1-1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmXYlm, m);
                FillBlock(l1, m1, l1+1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmXYlm, m);
                FillBlock(l1, m1, l1-1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmXYlm, m);
            }
        }

        m->AssembleBegin();
        m->AssembleEnd();
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
}
void YukawaPotential::FillMatrixGradY(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = basis.getOrder();
    int mmax = Ms.back();

    if (_isCentral) {
        // the same for each L-M block so cache it
        std::vector<complex> invR(N*N), invRR(N*N);
        BuildExpInvR(basis, N, _Z, _D, invR);
        for (auto& i : invR) i*=-_D;                 // from product rule
        BuildExpInvRR(basis, N, _Z, _D, invRR);

        // for each m-block (block rows)
        for (int m1 : Ms) {
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
                FillBlock(l1, m1, l1+1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmYYlm, m);
                FillBlock(l1, m1, l1-1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmYYlm, m);
                FillBlock(l1, m1, l1+1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmYYlm, m);
                FillBlock(l1, m1, l1-1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmYYlm, m);
            }
        }

        m->AssembleBegin();
        m->AssembleEnd();
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
}
void YukawaPotential::FillMatrixGradZ(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = basis.getOrder();
    int mmax = Ms.back();

    if (_isCentral) {
        // the same for each L-M block so cache it
        std::vector<complex> invR(N*N), invRR(N*N);
        BuildExpInvR(basis, N, _Z, _D, invR);
        for (auto& i : invR) i*=-_D;                 // from product rule
        BuildExpInvRR(basis, N, _Z, _D, invRR);


        for (int m1 : Ms) {
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) { 
                FillBlock(l1, m1, l1+1, m1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmZYlm, m);
                FillBlock(l1, m1, l1-1, m1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmZYlm, m);
            }
        }
        m->AssembleBegin();
        m->AssembleEnd();
    } else {
        assert(_isCentral && "only supporting central potentials.");
    } 
}




// utility functions
void BuildExpInvR(const Basis::BSpline& basis, int N, double Z, double D, std::vector<complex>& invR) {
    invR.resize(N*N);
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            invR[i + j*N] = invR[j + i*N] = basis.Integrate(i+1, j+1, [Z,D] (complex r) -> complex {
                return (-Z/r)*std::exp(-D*r);
            });
        }
    }
}
void BuildExpInvRR(const Basis::BSpline& basis, int N, double Z, double D, std::vector<complex>& invRR) {
    invRR.resize(N*N);
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            invRR[i + j*N] = invRR[j + i*N] = basis.Integrate(i+1, j+1, [Z,D] (complex r) -> complex {
                return (Z/r/r)*std::exp(-D*r);
            });
        }
    }
}

void FillBlock( int l1, int m1, int l2, int m2,
                int lmax, int mmax, int N, int order,
                const std::vector<int>& Ms,
                const std::vector<int>& mRows,
                const std::vector<complex>& invR,
                const std::vector<complex>& invRR,
                std::function<complex(int, int, int, int)> YlmYlm,
                Matrix m) {
    if (l2 <= lmax && l2 >= std::abs(m2) &&
        l2 >= 0 && l2 >= std::abs(m2) &&
        m2 <= mmax && m2 >= -mmax) {
        int m1Block = RowFrom(m1, Ms, mRows)/N;
        int m2Block = RowFrom(m2, Ms, mRows)/N;
        int blockRow = m1Block + (l1-std::abs(m1));
        int blockCol = m2Block + (l2-std::abs(m2));

        m->FillBandedBlock(order-1, N, blockRow, blockCol, 
        [=,&invR,&invRR](int row, int col) {
            int i = row % N, j = col % N;
            return (invR[i + j*N] + invRR[i + j*N])*YlmYlm(l1,m1,l2,m2);
        });
    }
}