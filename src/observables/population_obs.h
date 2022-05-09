#pragma once

#include "tdse/tdse.h"
#include "tdse/observable.h"
#include "maths/maths.h"

class PopulationObservable : public Observable {
    Matrix _S;
    Vector _psi;            // shortcut to wavefunction
    Vector _psi_temp;

    int _N, _lmax;
    int _eigen_state_nmax;
    std::vector<std::vector<Vector>> _states;       // (l,n)
    
    ASCII _file;
public:
    PopulationObservable(TDSE& tdse);

    void Startup();
    void Shutdown();
    void Compute(int it, double t, double  dt);
    void Flush();
};
