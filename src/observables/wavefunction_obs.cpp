
#include "observables/wavefunction_obs.h"

WavefunctionObservable::WavefunctionObservable(TDSE& tdse) : 
    Observable(tdse), _numGrid(0) {}

void WavefunctionObservable::SetNumGrid(int numGrid) {
    _numGrid = numGrid;
}

#include <iostream>
void WavefunctionObservable::Startup() {
    _psi = _tdse.Psi();
    _psi_grid = _MathLib.CreateVector(_numGrid);
    _grid.resize(_numGrid);

    double dx = (_tdse.Xmax() - _tdse.Xmin())/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        _grid[i] = _tdse.Xmin() + i*dx;

    _hdf5 = _MathLib.OpenHDF5(_output_filename, 'w');
}
void WavefunctionObservable::Shutdown() {
    _psi = nullptr;
    _psi_grid = nullptr;
    _hdf5 = nullptr;
}
void WavefunctionObservable::Compute(int it, double t, double dt) {
    std::stringstream ss;
    auto& basis = _tdse.Basis();
    _psi->Transform(_psi_grid, [=,&basis](const std::vector<complex>& coeff) {
        auto v = basis.FunctionEvaluate(_grid, coeff);
        return v;
    });
    
    _hdf5->PushGroup("vectors");
    for (int l = 0; l < _tdse.Lmax(); l++) {            // l-quantum number
        for (int i = 0; i < _numGrid; i++) {            // n-quantum number
            ss.str("");                                 // clear string stream
            ss << "psi[" << it << "]";

            _hdf5->WriteVector(ss.str().c_str(), _psi_grid);
            _hdf5->WriteAttribute(_psi_grid, "time", t);
        }
    }    
    _hdf5->PopGroup();
}

