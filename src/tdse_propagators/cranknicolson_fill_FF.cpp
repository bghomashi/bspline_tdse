
#include "tdse_propagators/cranknicolson.h"
#include "utility/index_manip.h"

void CrankNicolsonTDSE::FillFieldFree(Matrix& H0) {
    // int nl = _lmax + 2*_lmax*_mmax - _mmax*(_mmax+1); ?

    // TODO: decide on banded structure in non-central case
    // - for now assuming central
    Matrix temp = _MathLib.CreateMatrix(_dof, _dof, 2*_order-1);
    std::vector<complex> kinBlockStore(_N*_N);
    std::vector<complex> r2BlockStore(_N*_N);

    // Fill H0 with kinetic energy
    // derivative part is always the same (only depends on i,j)
    for (int i = 0; i < _N; i++) {
        for (int j = i; j < _N; j++) {
            kinBlockStore[i + j*_N] = kinBlockStore[j + i*_N] = _basis.Integrate(i+1, j+1, 1,1) /2.;
            r2BlockStore [i + j*_N] =  r2BlockStore[j + i*_N] = _basis.Integrate(i+1, j+1, [] (complex r) {
                return 1./r/r;
            });
        }
    }
    H0->FillBandedBlock(_order-1, _N, [=](int row, int col) {
        int i, j, l1, l2, m1, m2;
        ILMFrom(row, i, l1, m1, _N, _Ms, _mRows);
        ILMFrom(col, j, l2, m2, _N, _Ms, _mRows);
        // assert(l1 == l2 && m1 == m2);
        // kinetic energy and centrifugal term
        return kinBlockStore[i + j*_N] + 0.5*l1*(l1+1.)*r2BlockStore[i + j*_N];
    }); //, _N, _lmax, _mmax);

    for (auto& p : _potentials) {
        p->FillMatrix(_basis, temp, _N, _Ms, _mRows);    // get potential matrix
        _MathLib.AXPY(H0, 1., temp);                      // add potential to H0
    }
}