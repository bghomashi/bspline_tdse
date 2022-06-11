#include "tdse/tdse.h"
#include "observables/basis_obs.h"
#include <sstream>
#include <iomanip>
#include <limits.h>

BasisObservable::BasisObservable(TDSE& tdse) : Observable(tdse), _from(0), _to(INT_MAX) {}         

void BasisObservable::Startup(int it) {
    // dump the whole field into a text file
    if (_output_filename.length() > 0) {
        _txt_file = ASCII(_MathLib.OpenASCII(_output_filename, 'w'));
    }

    auto& basis = _tdse.Basis();
    int num_bsplines = basis.getNumBSplines();
    _to = std::min(num_bsplines-1, _to);
    std::vector<double> grid(_numGrid);
    std::vector<std::vector<complex>> splines(_to - _from + 1);

    double dx = (_tdse.Xmax() - _tdse.Xmin())/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        grid[i] = _tdse.Xmin() + i*dx;
        
    for (int bs = 0; bs < splines.size(); bs++)
        splines[bs] = basis.getBSpline(grid, bs + _from);

    // set up stringstream for formating and grab some shortcut to values we need
    std::stringstream ss;
    ss << std::setprecision(8) << std::scientific;
    
    // loop though the entire pulse
    for (int i = 0; i < _numGrid; i++) {
        ss.str("");
        ss << grid[i];
        for (auto& bs : splines) {
            ss << "\t" << std::real(bs[i]);
            ss << "\t" << std::imag(bs[i]);
        }
        ss << "\n";
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
void BasisObservable::SetFrom(int from) {
    _from = from;
} 
void BasisObservable::SetTo(int to) {
    _to = to;
}