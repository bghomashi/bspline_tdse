
#include "tdse_propagators/cranknicolson.h"
#include "utility/index_manip.h"



void CrankNicolsonTDSE::FillOverlap(Matrix& S) {
    std::vector<complex> overlapStore(_N*_N);
    for (int i = 0; i < _N; i++) {
        for (int j = i; j < _N; j++) {
            overlapStore[j + i*_N] = overlapStore[i + j*_N] = _basis.Integrate(i+1, j+1);
        }
    }
    S->FillBandedBlock(_order-1, _N, [=,&overlapStore](int row, int col) {
        int i, j, l1, l2, m1, m2;
        ILMFrom(row, i, l1, m1, _N, _Ms, _mRows);
        ILMFrom(col, j, l2, m2, _N, _Ms, _mRows);

        return overlapStore[i + j*_N];
    });
}
