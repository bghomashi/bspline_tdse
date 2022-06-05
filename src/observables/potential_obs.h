#pragma once

#include "tdse/observable.h"

class PotentialObservable : public Observable {
protected:
    ASCII _txt_file;
    int _nx;
public:
    PotentialObservable(TDSE& tdse);         
    
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it, double t, double dt);

    void SetNumGrid(int dx);
};