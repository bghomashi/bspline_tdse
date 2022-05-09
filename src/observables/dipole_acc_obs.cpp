#include "observables/dipole_acc_obs.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include "utility/logger.h"

#include "math_libs/petsc/petsc_lib.h"

DipoleAccObservable::DipoleAccObservable(TDSE& tdse) : Observable(tdse) {
}
void DipoleAccObservable::Startup() {
    // Build GradPotential Matrix
    auto& basis = _tdse.Basis();
    int Nmax = basis.getNumBSplines();
    int Lmax = _tdse.Lmax();
    int order = basis.getOrder();
    auto& Ms = _tdse.Ms();
    auto& mRows = _tdse.MRows();
    auto& potentials = _tdse.Potentials();
    auto polarization = _tdse.Polarization();
    
    _psi = _tdse.Psi();
    _psi_temp = _MathLib.CreateVector(_psi->Length());

    if (polarization[X]) {
        Matrix temp = _MathLib.CreateMatrix(_tdse.DOF(), _tdse.DOF(), 8*order-4);
        _gradPot[X] = _MathLib.CreateMatrix(_tdse.DOF(), _tdse.DOF(), 8*order-4);

        // Fill
        potentials[0]->FillMatrixGradX(basis, _gradPot[X], Nmax, Lmax, Ms, mRows);
        _gradPot[X]->Scale(-1.);
        for (int i = 1; i < potentials.size(); i++) {
            const auto& p = potentials[i]; 
            p->FillMatrixGradX(basis, temp, Nmax, Lmax, Ms, mRows);    // get potential matrix
            _MathLib.AXPY(_gradPot[X], -1., temp);                   // add potential to H0
        }
    }
    if (polarization[Y]) {
        Matrix temp = _MathLib.CreateMatrix(_tdse.DOF(), _tdse.DOF(), 8*order-4);
        _gradPot[Y] = _MathLib.CreateMatrix(_tdse.DOF(), _tdse.DOF(), 8*order-4);

        // Fill
        potentials[0]->FillMatrixGradY(basis, _gradPot[Y], Nmax, Lmax, Ms, mRows);
        _gradPot[Y]->Scale(-1.);
        for (int i = 1; i < potentials.size(); i++) {
            const auto& p = potentials[i]; 
            p->FillMatrixGradY(basis, temp, Nmax, Lmax, Ms, mRows);    // get potential matrix
            _MathLib.AXPY(_gradPot[Y], -1., temp);                   // add potential to H0
        }
    }
    if (polarization[Z]) {
        Matrix temp = _MathLib.CreateMatrix(_tdse.DOF(), _tdse.DOF(), 4*order-2);
        _gradPot[Z] = _MathLib.CreateMatrix(_tdse.DOF(), _tdse.DOF(), 4*order-2);

        // Fill
        potentials[0]->FillMatrixGradZ(basis, _gradPot[Z], Nmax, Lmax, Ms, mRows);
        _gradPot[Z]->Scale(-1.);
        for (int i = 1; i < potentials.size(); i++) {
            const auto& p = potentials[i]; 
            p->FillMatrixGradZ(basis, temp, Nmax, Lmax, Ms, mRows);    // get potential matrix
            _MathLib.AXPY(_gradPot[Z], -1., temp);                   // add potential to H0
        }
    }
    _txt_file = _MathLib.OpenASCII(_output_filename, 'w');
}

void DipoleAccObservable::Shutdown() {
    _psi = nullptr;
    _psi_temp = nullptr;
    _gradPot[X] = nullptr;
    _gradPot[Y] = nullptr;
    _gradPot[Z] = nullptr;
    _txt_file = nullptr;
}


void DipoleAccObservable::Compute(int it, double t, double dt) {
    complex dipole[DimIndex::NUM] = {0., 0., 0.};
    std::stringstream ss;

    if (_gradPot[X]) {
        _MathLib.Mult(_gradPot[X], _psi, _psi_temp);
        _MathLib.Dot(_psi, _psi_temp, dipole[X]);
    }
    if (_gradPot[Y]) {
        _MathLib.Mult(_gradPot[Y], _psi, _psi_temp);
        _MathLib.Dot(_psi, _psi_temp, dipole[Y]);
    }
    if (_gradPot[Z]) {
        _MathLib.Mult(_gradPot[Z], _psi, _psi_temp);
        _MathLib.Dot(_psi, _psi_temp, dipole[Z]);
    }

    // MatView(std::dynamic_pointer_cast<PetscMatrix>(_gradPotential)->_petsc_mat, 0);
    // exit(0);
    // VecView(std::dynamic_pointer_cast<PetscVector>(_psi)->_petsc_vec, 0);
    // exit(0);
    // VecView(std::dynamic_pointer_cast<PetscVector>(_psi_temp)->_petsc_vec, 0);
    //exit(0);
    ss << std::setprecision(8) << std::scientific;
    ss  << t << "\t" 
        << std::real(dipole[X]) << "\t"
        << std::real(dipole[Y]) << "\t"
        << std::real(dipole[Z]) << std::endl;
    _txt_file->Write(ss.str().c_str());
}

void DipoleAccObservable::Flush() {
    _txt_file->Flush();
}

int DipoleAccObservable::MemoryAlloced() const {
    int order = _tdse.Basis().getOrder();
    int mem = _tdse.DOF();
    if (_gradPot[X])
        mem += _tdse.DOF()*(8*order-4);
    if (_gradPot[Y])
        mem += _tdse.DOF()*(8*order-4);
    if (_gradPot[Z])
        mem += _tdse.DOF()*(4*order-2);
    
    return mem;
}