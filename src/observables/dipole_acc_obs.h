#pragma once

#include "tdse/tdse.h"
#include "tdse/observable.h"
#include "maths/maths.h"

class DipoleAccObservable : public Observable {
    Matrix _gradPot[DimIndex::NUM];
    Vector _psi;            // shortcut to wavefunction
    Vector _psi_temp;       // just from MatMult output

    ASCII _txt_file;
public:
    DipoleAccObservable(TDSE& tdse);

    int MemoryAlloced() const;
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it, double t, double dt);
};
