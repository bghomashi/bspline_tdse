#pragma once

#include "tdse/tdse.h"
#include <complex>
#include <map>

class CrankNicolsonTDSE : public TDSE {
    GMRESSolver _solver;

    Vector _psi_temp;
    Matrix _U0p, _U0m, _HI[DimIndex::NUM];
    Matrix _Up, _Um;
public:
    CrankNicolsonTDSE(MathLib& lib);
    void Initialize();
    bool DoStep(int it, double t, double dt);
    void Finish();

    void FillFieldFree(Matrix& m);
    void FillOverlap(Matrix& m);
    void FillInteractionX(Matrix& m);
    void FillInteractionY(Matrix& m);
    void FillInteractionZ(Matrix& m);
    void FillU0(Matrix& m);

    void DoCheckpoint();
    void DoObservables();
};
