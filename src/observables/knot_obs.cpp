#include "tdse/tdse.h"
#include "observables/knot_obs.h"
#include <sstream>
#include <iomanip>
#include <limits.h>

KnotObservable::KnotObservable(TDSE& tdse) : Observable(tdse) {}         

void KnotObservable::Startup(int it) {
    // dump the whole field into a text file
    if (_output_filename.length() > 0) {
        _txt_file = ASCII(_MathLib.OpenASCII(_output_filename, 'w'));
    }

    auto& basis = _tdse.Basis();
    auto& grid = basis.getGrid();
    
    // set up stringstream for formating and grab some shortcut to values we need
    std::stringstream ss;
    ss << std::setprecision(8) << std::scientific;
    
    // loop though the entire pulse
    for (int i = 0; i < grid.size(); i++) {
        ss.str("");
        ss << i << "\t"
           << grid[i] << "\n";
        _txt_file->Write(ss.str());
    }

    // done
    _txt_file = nullptr;
}



void KnotObservable::Flush() {}
void KnotObservable::Shutdown() {}
void KnotObservable::Compute(int it, double t, double dt) {}