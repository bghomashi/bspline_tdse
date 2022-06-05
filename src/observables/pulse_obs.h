#pragma once

#include "tdse/observable.h"

class PulseObservable : public Observable {
protected:
    ASCII _txt_file;
public:
    PulseObservable(TDSE& tdse);         
    
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it, double t, double dt);
};