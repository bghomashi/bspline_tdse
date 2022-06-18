#pragma once

#include "tdse/observable.h"

class KnotObservable : public Observable {
protected:
    ASCII _txt_file;
public:
    KnotObservable(TDSE& tdse);
        
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it, double t, double dt);
};