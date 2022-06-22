#pragma once

#include <string>
#include "tdse/tdse.h"
#include "tdse/observable.h"
#include "maths/maths.h"

// this while output the 3D density in cartesian coordinates. Use responsibly
class DensityObservable : public Observable {
    ASCII _txt_file;

    int _numGrid;
    std::vector<double> _grid;
public:
    DensityObservable(TDSE& tdse);

    void SetNumGrid(int numGrid);

    void Startup(int it);
    void Shutdown();
    void Compute(int it, double t, double dt);
};
