#include "tdse/tdse.h"
#include "observables/pulse_obs.h"
#include <sstream>
#include <iomanip>

PulseObservable::PulseObservable(TDSE& tdse) : Observable(tdse) {}         

void PulseObservable::Startup(int it) {
    // dump the whole field into a text file
    if (_output_filename.length() > 0) {
        _txt_file = ASCII(_MathLib.OpenASCII(_output_filename, 'w'));
    }

    // set up stringstream for formating and grab some shortcut to values we need
    std::stringstream ss;
    ss << std::setprecision(8) << std::scientific;
    double dt = _tdse.Timestep();
    double tmin = _tdse.Tmin();
    int NT = _tdse.NumTimeSteps();
    auto& field_x = _tdse.GetField(DimIndex::X);
    auto& field_y = _tdse.GetField(DimIndex::Y);
    auto& field_z = _tdse.GetField(DimIndex::Z);
    
    // loop though the entire pulse
    for (int i = 0; i < NT; i++) {

        ss.str("");
        ss << tmin + dt*i;
        if (_tdse.Polarization()[DimIndex::X])
            ss << "\t" << field_x[i];
        else 
            ss << "\t" << 0.0;
        if (_tdse.Polarization()[DimIndex::Y])
            ss << "\t" << field_y[i];
        else 
            ss << "\t" << 0.0;
        if (_tdse.Polarization()[DimIndex::Z])
            ss << "\t" << field_z[i];
        else 
            ss << "\t" << 0.0;
        ss << "\n";
        _txt_file->Write(ss.str());
    }

    // done
    _txt_file = nullptr;
}
void PulseObservable::Shutdown() {}
void PulseObservable::Compute(int it, double t, double dt) {}
void PulseObservable::Flush() {}