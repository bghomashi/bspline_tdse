#include "observables/density_observable.h"
#include "core/utility/logger.h"
#include "core/utility/index_manip.h"
#include "core/utility/spherical_harmonics.h"
#include <iomanip>
#include <cmath>

DensityObservable::DensityObservable(TDSE& tdse) : 
    Observable(tdse), _numGrid(0) {}

void DensityObservable::SetNumGrid(int numGrid) {
    _numGrid = numGrid;
}

void DensityObservable::Startup(int start_it) {
    // same grid in all dimensions
    _grid.resize(_numGrid);

    // make the grid symmetric [-max, max]
    double dx = 2.*_tdse.Xmax()/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        _grid[i] = -_tdse.Xmax() + i*dx;

}
void DensityObservable::Shutdown() {
    _txt_file = nullptr;
}
void DensityObservable::Compute(int it, double t, double dt) {
    _txt_file = _MathLib.OpenASCII("wavefunction_"+std::to_string(it)+".txt", 'w');

    std::vector<complex> psi(_tdse.DOF()), n_block_coeff;
    std::vector<complex> temp(_grid.size());
    std::vector<complex> out(_grid.size()*_grid.size()*_grid.size(), 0);
    std::stringstream ss;
    auto& basis = _tdse.Basis();
    int N = basis.getNumBSplines();
    int lmax = _tdse.Lmax();
    auto& Ms = _tdse.Ms();
    auto& MRows = _tdse.MRows();
    

    _tdse.Psi()->CopyTo(psi);                                                                   // copy entire wf into vector

    for (int m : Ms) {                                                                              // for each m
        for (int l = std::abs(m); l <= lmax; l++) {                                                 // for each l
            int start = RowFrom(0, m, l, N, Ms, MRows);                                             // first index of this block
            n_block_coeff = std::vector<complex>(psi.begin() + start, psi.begin() + (start+N)); // copy n_block coeffs
            
            for (int i = 0; i < _numGrid; i++) {
                for (int j = 0; j < _numGrid; j++) {
                    for (int k = 0; k < _numGrid; k++) {
                        double x = _grid[i];
                        double y = _grid[j];
                        double z = _grid[k];

                        double r = std::sqrt(x*x + y*y + z*z);
                        double theta = (x == 0 && y == 0 && z == 0 ? 0.0 : std::atan2(std::sqrt(x*x + y*y), z));
                        double phi = (x == 0 && y == 0 ? 0.0 : std::atan2(y, x));

                        out[i + j*N + k*N*N] += basis.FunctionEvaluate(r, n_block_coeff)*Ylm(l, m, theta, phi); // evaluate on grid this chuck on the grid
                    }
                }
            }            
        }  
    }       

    for (int i = 0; i < _numGrid; i++) {
        for (int j = 0; j < _numGrid; j++) {
            for (int k = 0; k < _numGrid; k++) {  
                ss.str("");
                ss << _grid[i] << "\t"; // x
                ss << _grid[j] << "\t"; // y
                ss << _grid[k] << "\t"; // z
                ss << std::real(out[i]) << "\t";
                ss << std::imag(out[i]) << "\t";
                ss << std::abs(out[i])*std::abs(out[i]) << "\n";
                _txt_file->Write(ss.str());  
            }
        }
    }

    _txt_file = nullptr;
}

