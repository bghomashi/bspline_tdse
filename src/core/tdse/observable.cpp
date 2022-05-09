#include "tdse/observable.h"
#include "tdse/tdse.h"

#include <iostream>

Observable::Observable(TDSE& tdse) : _tdse(tdse), _MathLib(tdse.MathLibrary()) {
    _compute_period_in_iterations = 1;
}
void Observable::DoObservable(int it, double t, double dt) {
    if (it % _compute_period_in_iterations == 0) {
        Compute(it, t, dt);
    }
}

void Observable::SetComputePeriod(int iterations) {
    _compute_period_in_iterations = iterations;
}
void Observable::SetFilename(const std::string& filename) {
    _output_filename = filename;
}

int Observable::MemoryAlloced() const {
    return 0;
}