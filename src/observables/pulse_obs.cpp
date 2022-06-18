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
    auto& pulses = _tdse.Pulses();
    // auto& A_field_x = _tdse.GetField(DimIndex::X);
    // auto& A_field_y = _tdse.GetField(DimIndex::Y);
    // auto& A_field_z = _tdse.GetField(DimIndex::Z);
    
    std::vector<Vec3> A_field(NT, {0,0,0});
    std::vector<Vec3> E_field(NT, {0,0,0});

    for (auto& pulse : pulses) {
        for (int i = 0; i < NT; i++) {
            A_field[i] += pulse->A(tmin + dt*i);
            E_field[i] += pulse->E(tmin + dt*i);
        }
    }

    // loop though the entire pulse
    for (int i = 0; i < NT; i++) {

        ss.str("");
        ss << tmin + dt*i;
        // A-field
        /*
        if (_tdse.Polarization()[DimIndex::X])
            ss << "\t" << A_field[i].x;
        else 
            ss << "\t" << 0.0;
        if (_tdse.Polarization()[DimIndex::Y])
            ss << "\t" << A_field[i].y;
        else 
            ss << "\t" << 0.0;
        if (_tdse.Polarization()[DimIndex::Z])
            ss << "\t" << A_field[i].z;
        else 
            ss << "\t" << 0.0;
        */

        ss << "\t" << A_field[i].x;
        ss << "\t" << A_field[i].y;
        ss << "\t" << A_field[i].z;
        ss << "\t" << E_field[i].x;
        ss << "\t" << E_field[i].y;
        ss << "\t" << E_field[i].z;

        ss << "\n";
        _txt_file->Write(ss.str());
    }

    // done
    _txt_file = nullptr;
}
void PulseObservable::Shutdown() {}
void PulseObservable::Compute(int it, double t, double dt) {}
void PulseObservable::Flush() {}