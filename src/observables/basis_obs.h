#pragma once

#include "tdse/observable.h"

class BasisObservable : public Observable {
protected:
    ASCII _txt_file;
    int _numGrid, _to, _from;
    std::vector<double> _grid;
public:
    BasisObservable(TDSE& tdse);
    
    void SetFrom(int from);  
    void SetTo(int to);  
    void SetNumGrid(int numGrid);  
    
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it, double t, double dt);
};