#pragma once

#include "maths/maths.h"
#include <memory>

// Observable-objects live inside the TDSE-object
// They have access to the other components of the TDSE-class
class TDSE;

class Observable {
protected:
    TDSE& _tdse;                        // parent TDSE-object
    MathLib& _MathLib;                  // - shortcut

    int _compute_period_in_iterations;
    std::string _output_filename;
public:
    typedef std::shared_ptr<Observable> Ptr_t;

    // since the observables will be aggregated inside 
    // the TDSE object they will be destroyed with (before) it.
    // So we do not need to worry about the TDSE& reference becoming invalid.
    Observable(TDSE& tdse);         

    void SetComputePeriod(int iterations);
    void DoObservable(int it, double t, double dt);
    void SetFilename(const std::string& filename);
    
    virtual void Flush() {};
    virtual int MemoryAlloced() const;
    virtual void Startup(int it) = 0;
    virtual void Shutdown() = 0;
    virtual void Compute(int it, double t, double dt) = 0;
};