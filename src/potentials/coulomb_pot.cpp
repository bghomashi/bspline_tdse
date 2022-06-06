#include "potentials/coulomb_pot.h"
#include "utility/index_manip.h"
#include "utility/logger.h"
#include "utility/spherical_harmonics.h"

#include <iostream>
#include <functional>

static void BuildInvR(const Basis::BSpline& basis, int N, double Z, std::vector<complex>& invR);
static void BuildInvRR(const Basis::BSpline& basis, int N, double Z, std::vector<complex>& invRR);
static void FillBlock( int l1, int m1, int l2, int m2,
                int lmax, int mmax, int N, int order,
                const std::vector<int>& Ms,
                const std::vector<int>& mRows,
                const std::vector<complex>& invRR,
                std::function<complex(int, int, int, int)> YlmYlm,
                Matrix m);



CoulombPotential::CoulombPotential() : _Z(0), _x(0), _y(0), _z(0) {}
CoulombPotential::CoulombPotential(double Z) : _Z(Z), _x(0), _y(0), _z(0) {}
CoulombPotential::CoulombPotential(double Z, double x, double y, double z) :
    _Z(Z), _x(x), _y(y), _z(z) {
    if (z != 0) {
        _isCentral = false;
        if (x != 0 || y != 0)
            _isAxial = false;
    }
}

double CoulombPotential::operator() (double x, double y, double z) const {
    double r = std::sqrt(x*x + y*y + z*z);
    return 1./r;
}
// FIX: maybe this way of thinking about it is not very efficient
void CoulombPotential::FillMatrix(const Basis::BSpline& basis, Matrix m, int N, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = basis.getOrder();

    if (_isCentral) {
        // the same for each L-M block so cache it
        std::vector<complex> invR(N*N);
        BuildInvR(basis, N, _Z, invR);

        m->FillBandedBlock(order-1, N, [&](int row, int col) {
            int i, j, l1, l2, m1, m2;
            ILMFrom(row, i, l1, m1, N, Ms, mRows);
            ILMFrom(col, j, l2, m2, N, Ms, mRows);
            return invR[i + j*N];
        }); //, N, lmax, mmax);
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
        
    // if axial:
    // for each total_row r:
    //     for each total_col c:
    //         int i = r % N;
    //         int j = c % N;
    //         int  = r / N;
    //         for each l:
    //             basis.integrate(i,j, 0, R,  
    //                 [](complex r) -> complex {
    //                     return r^l/R^l+1;
    //                 }) + 
    //             basis.integrate(i,j, R, rmax,
    //                 [](complex r) -> complex {
    //                     return R^l/r^l+1;
    //                 }) *
    //             TripleIntegral(l1,l,l2, -m1,0,m2);
}

void CoulombPotential::FillMatrixGradX(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = basis.getOrder();
    // if we get here this is a full 3D calculations
    // so Ms = [-mmax,mmax]
    int mmax = Ms.back();

    if (_isCentral) {
        // the same for each L-M block so cache it
        std::vector<complex> invRR(N*N);
        BuildInvRR(basis, N, _Z, invRR);

        for (int m1 : Ms) {
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
                FillBlock(l1, m1, l1+1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmXYlm, m);
                FillBlock(l1, m1, l1-1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmXYlm, m);
                FillBlock(l1, m1, l1+1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmXYlm, m);
                FillBlock(l1, m1, l1-1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmXYlm, m);
            }
            // int m1Block = mRows[i]/N;                     // number of l-blocks to skip (rows)
            // int m2;
            // int m2Block;                                    // number of l-blocks to skip (cols)
            // int blockRow, blockCol, l2;

            // if (m1+1 <= mmax) {                            // then there is a coupling
            //     m2 = m1+1;
            //     m2Block = RowFrom(m2, Ms, mRows)/N;
            //     // for each l-block (block rows)
            //     for (int l1 = std::abs(m1); l1 <= lmax; l1++) {         // m1=m2-1, l1=l2-1
            //         // check one l-block up

            //         if (l1+1 <= lmax && l1+1 >= std::abs(m2)) {
            //             l2 = l1+1;
            //             blockRow = m1Block + (l1-std::abs(m1));
            //             blockCol = m2Block + (l2-std::abs(m2));

            //             m->FillBandedBlock(order-1, N, blockRow, blockCol, 
            //             [=,&invRR](int row, int col) {
            //                 int i = row % N, j = col % N;
            //                 return invRR[i + j*N]*YlmXYlm(l1,m1,l2,m2);
            //             });
            //         }

            //         // check one l-block down
            //         if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {            // m1=m2-1, l1=l2+1
            //             l2 = l1-1;
            //             blockRow = m1Block + (l1-std::abs(m1));
            //             blockCol = m2Block + (l2-std::abs(m2));

            //             m->FillBandedBlock(order-1, N, blockRow, blockCol, 
            //             [=,&invRR](int row, int col) {
            //                 int i = row % N, j = col % N;
            //                 return invRR[i + j*N]*YlmXYlm(l1,m1,l2,m2);
            //             });
            //         }
            //     }
            // }
            // if (m1-1 >= -mmax) {                           // then there is a coupling
            //     m2 = m1-1;
            //     m2Block = RowFrom(m2, Ms, mRows)/N;
            //     // for each l-block (block rows)
            //     for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
            //         // check one l-block up
            //         if (l1+1 <= lmax && l1+1 >= std::abs(m2)) {         // m1=m2+1, l1=l2-1
            //             l2 = l1+1;
            //             blockRow = m1Block + (l1-std::abs(m1));
            //             blockCol = m2Block + (l2-std::abs(m2));

            //             m->FillBandedBlock(order-1, N, blockRow, blockCol, 
            //             [=,&invRR](int row, int col) {
            //                 int i = row % N, j = col % N;
            //                 return invRR[i + j*N]*YlmXYlm(l1,m1,l2,m2);
            //             });
            //         }
            //         // check one l-block down
            //         if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {            // m1=m2+1, l1=l2+1
            //             l2 = l1-1;
            //             blockRow = m1Block + (l1-std::abs(m1));
            //             blockCol = m2Block + (l2-std::abs(m2));

            //             m->FillBandedBlock(order-1, N, blockRow, blockCol, 
            //             [=,&invRR](int row, int col) {
            //                 int i = row % N, j = col % N;
            //                 return invRR[i + j*N]*YlmXYlm(l1,m1,l2,m2);
            //             });
            //         }
            //     }
            
        }

        m->AssembleBegin();
        m->AssembleEnd();
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
}

void CoulombPotential::FillMatrixGradY(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = basis.getOrder();
    int mmax = Ms.back();

    if (_isCentral) {
        // the same for each L-M block so cache it
        std::vector<complex> invRR(N*N);
        BuildInvRR(basis, N, _Z, invRR);

        // for each m-block (block rows)
        for (int m1 : Ms) {
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
                FillBlock(l1, m1, l1+1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmYYlm, m);
                FillBlock(l1, m1, l1-1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmYYlm, m);
                FillBlock(l1, m1, l1+1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmYYlm, m);
                FillBlock(l1, m1, l1-1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmYYlm, m);
            }
        }

        m->AssembleBegin();
        m->AssembleEnd();
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
}

void CoulombPotential::FillMatrixGradZ(const Basis::BSpline& basis, Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = basis.getOrder();
    int mmax = Ms.back();

    if (_isCentral) {
        // the same for each L-M block so cache it
        std::vector<complex> invRR(N*N);
        BuildInvRR(basis, N, _Z, invRR);
        for (int m1 : Ms) {
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) { 
                FillBlock(l1, m1, l1+1, m1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmZYlm, m);
                FillBlock(l1, m1, l1-1, m1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmZYlm, m);
            }
        }
        m->AssembleBegin();
        m->AssembleEnd();
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
}



// utility functions
void BuildInvR(const Basis::BSpline& basis, int N, double Z, std::vector<complex>& invR) {
    invR.resize(N*N);
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            invR[i + j*N] = invR[j + i*N] = basis.Integrate(i+1, j+1, [Z] (complex r) -> complex {
                return -Z/r;
            });
        }
    }
}
void BuildInvRR(const Basis::BSpline& basis, int N, double Z, 
    std::vector<complex>& invRR) {
    invRR.resize(N*N);
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            invRR[i + j*N] = invRR[j + i*N] = basis.Integrate(i+1, j+1, [Z] (complex r) -> complex {
                return Z/r/r;
            });
        }
    }
}

void FillBlock( int l1, int m1, int l2, int m2,
                int lmax, int mmax, int N, int order,
                const std::vector<int>& Ms,
                const std::vector<int>& mRows,
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
        [=,&invRR](int row, int col) {
            int i = row % N, j = col % N;
            return invRR[i + j*N]*YlmYlm(l1,m1,l2,m2);
        });
    }
}