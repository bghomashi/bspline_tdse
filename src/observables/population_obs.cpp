#include "observables/population_obs.h"
#include "utility/logger.h"
#include "utility/index_manip.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include "math_libs/petsc/petsc_lib.h"

PopulationObservable::PopulationObservable(TDSE& tdse) : Observable(tdse) {
}
void PopulationObservable::Startup(int start_it) {
    auto& basis = _tdse.Basis();
    int order = basis.getOrder();
    _N = basis.getNumBSplines();
    _lmax = _tdse.Lmax();
    std::string eigen_state_filename = _tdse.GetInitialStateFile();
    _eigen_state_nmax = _tdse.GetInitialStateNmax();
    _eigen_state_lmax = _tdse.GetInitialStateLmax();

    _psi = _tdse.Psi();
    _S = _MathLib.CreateMatrix(_N, _N, 2*order-1);

    _psi_temp = _MathLib.CreateVector(_N);
   
    // storage for eigenstates
    // _states.resize(_eigen_state_nmax);
    // for (int l = 0; l < _eigen_state_nmax; l++)
    //     _states[l].resize(_eigen_state_nmax - l);
    
    // for (int l = 0; l < _eigen_state_nmax; l++)
    //     for (int n = 0; n < _eigen_state_nmax - l; n++)
    //         _states[l][n] = _MathLib.CreateVector(_N);


    auto hdf5 = _MathLib.OpenHDF5(eigen_state_filename, 'r');
    hdf5->PushGroup("vectors");

    // std::stringstream name_ss;
    // for (int l = 0; l < _eigen_state_nmax; l++) {
    //     for (int n = l+1; n <= _eigen_state_lmax; n++) {
    //         name_ss.str("");                                         // clear string stream
    //         name_ss << "(" << n << ", " << l << ")";     // name of state
    //         hdf5->ReadVector(name_ss.str().c_str(), _states[l][n-l-1]);     
    //     }
    // }
    
    hdf5->PopGroup();

    // Fill overlap matrix
    std::vector<complex> overlapStore(_N*_N);
    for (int i = 0; i < _N; i++) {
        for (int j = i; j < _N; j++) {
            overlapStore[j + i*_N] = overlapStore[i + j*_N] = basis.Integrate(i+1, j+1);
        }
    }
    _S->FillBandedBlock(order-1, _N, [=,&overlapStore](int row, int col) {
        int i, j;
        i = row % _N;  
        j = col % _N;  

        return overlapStore[i + j*_N];
    });




    if (_output_filename.length() > 0) {
        _file = ASCII(_MathLib.OpenASCII(_output_filename, 'w'));
    }
}
void PopulationObservable::Shutdown() {
    complex pop;
    std::stringstream ss;
    std::string eigen_state_filename = _tdse.GetInitialStateFile();

    ss << std::setprecision(8) << std::scientific;
    int min_lmax = std::min(_lmax, _eigen_state_lmax);

    int i = 0;
    std::vector<int> Ms = _tdse.Ms();
    std::vector<int> mRows = _tdse.MRows();


    auto hdf5 = _MathLib.OpenHDF5(eigen_state_filename, 'r');
    hdf5->PushGroup("vectors");


    Vector eigen_state = _MathLib.CreateVector(_N);

    if (_file)
        _file->Write("(n, l, m)\n");
    else
        Log::info("\n");
    for (auto m : Ms) {
        if (std::abs(m) >= _eigen_state_lmax)   // we do not have this state to project on to
            continue;                       // so skip
        for (int l = std::abs(m); l <= min_lmax; l++) {
            // get a subvector
            int start = RowFrom(m, Ms, mRows) + (l-std::abs(m))*_N;
            Vector ml_block = _psi->GetSubVector(start, start+_N);
            
            for (int n = l+1; n <= _eigen_state_nmax; n++) {
                // load vector from state file
                ss.str("");                                         // clear string stream
                ss << "(" << n << ", " << l << ")";     // name of state
                hdf5->ReadVector(ss.str().c_str(), eigen_state);    

                // project
                _MathLib.Mult(_S, eigen_state, _psi_temp);
                _MathLib.Dot(ml_block, _psi_temp, pop);
                
                ss.str("");
                ss << "(" << n << ", " << l << ", " << m << "): ";
                ss << std::abs(pop)*std::abs(pop) << "\t" 
                   << std::real(pop) << "\t" 
                   << std::imag(pop) << "\n";

                if (_file)
                    _file->Write(ss.str().c_str());
                else
                    Log::info(ss.str());
            }
            _psi->RestoreSubVector(ml_block);

            if (_file)
                _file->Write("\n");
            else
                Log::info("\n");
        }
        if (_file)
            _file->Write("\n");
        else
            Log::info("\n");
    }

    hdf5->PopGroup();

    // for (auto& l : _states) {
    //     for (auto& n : l)
    //         n = 0;
    // }

    hdf5 = nullptr;
    eigen_state = nullptr;
    _psi = nullptr;
    _psi_temp = nullptr;
    _S = nullptr;
    _file = nullptr;
}
void PopulationObservable::Compute(int it, double t, double dt) {
}


void PopulationObservable::Flush() {
    _file->Flush();
}