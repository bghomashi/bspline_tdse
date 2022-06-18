#pragma once

#include <string>
#include "tdse/tdse.h"
#include "tdse/observable.h"
#include "maths/maths.h"

class DebugEigenstatesObservable : public Observable {
    ASCII _txt_file;

    int _numGrid, _nmax;
    std::vector<double> _grid;
public:
    DebugEigenstatesObservable(TDSE& tdse);

    void SetNumGrid(int numGrid);
    void SetNMax(int nmax);
    void Startup(int it);
    void Shutdown();
    void Compute(int it, double t, double dt);
};
