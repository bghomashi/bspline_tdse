#pragma once

#include "tdse/tdse.h"
#include "tdse/observable.h"
#include "maths/maths.h"

class NormObservable : public Observable {
    Matrix _S;
    Vector _psi;            // shortcut to wavefunction
    Vector _psi_temp;

    ASCII _file;
public:
    NormObservable(TDSE& tdse);

    void Startup();
    void Shutdown();
    void Compute(int it, double t, double  dt);
    void Flush();
};
