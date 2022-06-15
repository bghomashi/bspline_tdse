#pragma once

#include <string>
#include "tdse/tdse.h"
#include "tdse/observable.h"
#include "maths/maths.h"

class DebugWavefunctionObservable : public Observable {
    ASCII _txt_file;

    int _numGrid;
    std::vector<double> _grid;
public:
    DebugWavefunctionObservable(TDSE& tdse);

    void SetNumGrid(int numGrid);

    void Startup(int it);
    void Shutdown();
    void Compute(int it, double t, double dt);
};
