#include "observables/debug_wavefunction_obs.h"
#include "core/utility/logger.h"
#include <iomanip>

DebugWavefunctionObservable::DebugWavefunctionObservable(TDSE& tdse) : 
    Observable(tdse), _numGrid(0) {}

void DebugWavefunctionObservable::SetNumGrid(int numGrid) {
    _numGrid = numGrid;
}

#include <iostream>
void DebugWavefunctionObservable::Startup(int start_it) {
    _grid.resize(_numGrid);

    double dx = (_tdse.Xmax() - _tdse.Xmin())/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        _grid[i] = _tdse.Xmin() + i*dx;

}
void DebugWavefunctionObservable::Shutdown() {
    _txt_file = _MathLib.OpenASCII("wavefunction_final.txt", 'w');

    std::vector<complex> psi(_tdse.DOF()), n_block_coeff;
    std::vector<complex> out(_grid.size(), 0), temp(_grid.size());
    std::stringstream ss;
    auto& basis = _tdse.Basis();
    
    _tdse.Psi()->CopyTo(psi);

    int N = basis.getNumBSplines();

    for (int i = 0; i < _tdse.DOF(); i += N) {
        n_block_coeff = std::vector<complex>(psi.begin() + i, psi.begin() + (i+N));             // copy n_block coeffs
        temp = basis.FunctionEvaluate(_grid, n_block_coeff);                                        // evaluate on grid
        
        for (int i = 0; i < _grid.size(); i++)                                                      // sum the square at all the grid points                  
            out[i] += std::abs(temp[i])*std::abs(temp[i]);
    }
    ss << std::setprecision(8) << std::scientific;
    for (int i = 0; i < _numGrid; i++) {            // n-quantum number
        ss.str("");
        ss << _grid[i] << "\t" 
           << std::real(out[i]) << "\n";
        _txt_file->Write(ss.str());
    }


    _txt_file = nullptr;
}
void DebugWavefunctionObservable::Compute(int it, double t, double dt) {
    _txt_file = _MathLib.OpenASCII("wavefunction_"+std::to_string(it)+".txt", 'w');

    std::vector<complex> psi(_tdse.DOF()), n_block_coeff;
    std::vector<complex> out(_grid.size(), 0), temp(_grid.size());
    std::stringstream ss;
    auto& basis = _tdse.Basis();
    
    _tdse.Psi()->CopyTo(psi);

    int N = basis.getNumBSplines();

    for (int i = 0; i < _tdse.DOF(); i += N) {
        n_block_coeff = std::vector<complex>(psi.begin() + i, psi.begin() + (i+N));             // copy n_block coeffs
        temp = basis.FunctionEvaluate(_grid, n_block_coeff);                                        // evaluate on grid
        
        for (int i = 0; i < _grid.size(); i++)                                                      // sum the square at all the grid points                  
            out[i] += std::abs(temp[i])*std::abs(temp[i]);
    }

    for (int i = 0; i < _numGrid; i++) {            // n-quantum number
        ss.str("");
        ss << _grid[i] << "\t" 
           << std::real(out[i]) << "\n";
        _txt_file->Write(ss.str());
    }


    _txt_file = nullptr;
}

