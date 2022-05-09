#pragma once

#include <string>
#include "tdse/tdse.h"
#include "tdse/observable.h"
#include "maths/maths.h"

class WavefunctionObservable : public Observable {
    Vector _psi;            // shortcut to wavefunction
    Vector _psi_grid;       // holds wavefunction in real-space 
    HDF5 _hdf5;

    int _numGrid;
    std::vector<double> _grid;
public:
    WavefunctionObservable(TDSE& tdse);

    void SetNumGrid(int numGrid);

    void Startup();
    void Shutdown();
    void Compute(int it, double t, double dt);
};
