#include "tdse/tdse.h"
#include "observables/potential_obs.h"
#include <sstream>
#include <iomanip>

PotentialObservable::PotentialObservable(TDSE& tdse) : Observable(tdse) {}         

void PotentialObservable::Startup(int it) {
    // dump the whole field into a text file
    if (_output_filename.length() > 0) {
        _txt_file = ASCII(_MathLib.OpenASCII(_output_filename, 'w'));
    }

    // set up stringstream for formating and grab some shortcut to values we need
    std::stringstream ss;
    ss << std::setprecision(8) << std::scientific;
    auto& potentials = _tdse.Potentials();
    double xmin = _tdse.Xmin();
    double xmax = _tdse.Xmax();
    double dx = (xmax - xmin)/(_nx - 1.);
    
    // this assumes all potentials are central!
    // loop though the entire pulse
    for (int i = 0; i < _nx; i++) {
        double r = xmin + dx*i;
        double p = 0;
        for (auto& pot : potentials)
            p += (*pot)(r,0,0);

        ss.str("");
        ss << r << "\t" << p << "\n";
        _txt_file->Write(ss.str());
    }

    // done
    _txt_file = nullptr;
}


void PotentialObservable::SetNumGrid(int nx) {
    _nx = nx;
}

void PotentialObservable::Shutdown() {}
void PotentialObservable::Compute(int it, double t, double dt) {}
void PotentialObservable::Flush() {
}