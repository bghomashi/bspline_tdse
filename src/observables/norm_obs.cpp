#include "observables/norm_obs.h"
#include "utility/logger.h"
#include "utility/file_exists.h"
#include <iostream>
#include <sstream>

NormObservable::NormObservable(TDSE& tdse) : Observable(tdse) {
}
void NormObservable::Startup(int start_it) {
    auto& basis = _tdse.Basis();
    int order = basis.getOrder();
    int NBsplines = basis.getNumBSplines();

    _S = _MathLib.CreateMatrix(_tdse.DOF(), _tdse.DOF(), 2*order-1);
    _psi = _tdse.Psi();
    _psi_temp = _MathLib.CreateVector(_psi->Length());

    // Fill overlap matrix
    std::vector<complex> overlapStore(NBsplines*NBsplines);
    for (int i = 0; i < NBsplines; i++) {
        for (int j = i; j < NBsplines; j++) {
            overlapStore[j + i*NBsplines] = overlapStore[i + j*NBsplines] = basis.Integrate(i+1, j+1);
        }
    }
    _S->FillBandedBlock(order-1, NBsplines, [=,&overlapStore](int row, int col) {
        int i, j;
        i = row % NBsplines;  
        j = col % NBsplines;  

        return overlapStore[i + j*NBsplines];
    });

    if (start_it > 0 && file_exists(_output_filename)) {
        std::vector<complex> norm((start_it+1)/_compute_period_in_iterations);
        std::stringstream line;
        
        // first read all the values back that we want to keep
        _file = _MathLib.OpenASCII(_output_filename, 'r');
        for (int i = 0; i < norm.size(); i++) {
            line.str(_file->ReadLine());
            line >> norm[i];
        }

        // now clear the file and write them all back
        _file = _MathLib.OpenASCII(_output_filename, 'w');
        for (int i = 0; i < norm.size(); i++) {
            line.str("");
            line << norm[i] << std::endl;
            _file->Write(line.str().c_str());
        }
    } else {
        _file = _MathLib.OpenASCII(_output_filename, 'w');
    }
}
void NormObservable::Shutdown() {
    _psi = nullptr;
    _psi_temp = nullptr;
    _S = nullptr;
    _file = nullptr;
}
void NormObservable::Compute(int it, double t, double dt) {
    complex norm;
    _MathLib.Mult(_S, _psi, _psi_temp);
    _MathLib.Dot(_psi, _psi_temp, norm);
    if (_file) {
        std::stringstream ss;
        ss << norm << "\n";
        _file->Write(ss.str().c_str());
    } else {
        std::stringstream ss;
        ss << "Norm: " << norm << "\n";
        Log::info(ss.str());
    }
}


void NormObservable::Flush() {
    _file->Flush();
}