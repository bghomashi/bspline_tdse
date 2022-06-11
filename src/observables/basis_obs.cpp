#include "tdse/tdse.h"
#include "observables/basis_obs.h"
#include <sstream>
#include <iomanip>

BasisObservable::BasisObservable(TDSE& tdse) : Observable(tdse) {}         

void BasisObservable::Startup(int it) {
    // dump the whole field into a text file
    if (_output_filename.length() > 0) {
        _txt_file = ASCII(_MathLib.OpenASCII(_output_filename, 'w'));
    }

    std::vector<double> grid(_numGrid);
    std::vector<complex> spline(_numGrid);

    double dx = (_tdse.Xmax() - _tdse.Xmin())/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        grid[i] = _tdse.Xmin() + i*dx;

    auto& basis = _tdse.Basis();
    int num_bsplines = basis.getNumBSplines();
    // set up stringstream for formating and grab some shortcut to values we need
    std::stringstream ss;
    ss << std::setprecision(8) << std::scientific;
    
    // loop though the entire pulse
    for (int i = 0; i < num_bsplines; i++) {
        spline = basis.getBSpline(grid, i);

        ss.str("");
        ss << grid[i] << "\t";
        ss << std::real(spline[i]) << "\t";
        ss << std::real(spline[i]) << "\n";
        
        _txt_file->Write(ss.str());
    }

    // done
    _txt_file = nullptr;
}
void BasisObservable::Shutdown() {}
void BasisObservable::Compute(int it, double t, double dt) {}
void BasisObservable::Flush() {}
void BasisObservable::SetNumGrid(int numGrid) {
    _numGrid = numGrid;
}