#include "observables/debug_eigenstates.h"
#include "core/utility/logger.h"
#include "core/tdse/tdse.h"
#include <iomanip>

DebugEigenstatesObservable::DebugEigenstatesObservable(TDSE& tdse) : 
    Observable(tdse), _numGrid(0) {}

void DebugEigenstatesObservable::SetNumGrid(int numGrid) {
    _numGrid = numGrid;
}

void DebugEigenstatesObservable::Startup(int start_it) {
    _grid.resize(_numGrid);

    double dx = (_tdse.Xmax() - _tdse.Xmin())/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        _grid[i] = _tdse.Xmin() + i*dx;

    std::stringstream name_ss, ss;
    std::vector<complex> temp(_grid.size()), out(_grid.size());
    auto& basis = _tdse.Basis();
    int N = basis.getNumBSplines();
    Vector eigen = _MathLib.CreateVector(N);        // vector from hdf5 file

    auto hdf5 = _MathLib.OpenHDF5(_tdse.GetInitialStateFile(), 'r');
    hdf5->PushGroup("vectors");

    ss << std::setprecision(8) << std::scientific;
    for (int n = 1; n <= _nmax; n++) {
        for (int l = 0; l < n; l++) {
            name_ss.str("");                                         // clear string stream
            name_ss << "(" << n << ", " << l << ")";     // name of state

            _txt_file = _MathLib.OpenASCII("eigenstate_"+name_ss.str()+".txt", 'w');

            hdf5->ReadVector(name_ss.str().c_str(), eigen);           // read in the state
            eigen->CopyTo(temp);
            out = basis.FunctionEvaluate(_grid, temp);  

            for (int i = 0; i < _numGrid; i++) {            // n-quantum number
                ss.str("");
                ss << _grid[i] << "\t" <<
                      std::real(out[i]) << "\n";
                _txt_file->Write(ss.str());
            }                                      // evaluate on grid
        }
    }
    hdf5->PopGroup();

    hdf5 = nullptr;
    eigen = nullptr;
    _txt_file = nullptr;
}

void DebugEigenstatesObservable::SetNMax(int nmax) {
    _nmax = nmax;
}
void DebugEigenstatesObservable::Compute(int it, double t, double dt) {}
void DebugEigenstatesObservable::Shutdown() {}