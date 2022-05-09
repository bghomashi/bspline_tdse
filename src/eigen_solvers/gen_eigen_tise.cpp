
#include "eigen_solvers/gen_eigen_tise.h"
#include "utility/logger.h"

#include <string>
#include <sstream>

GeneralizeEigenvalueTISE::GeneralizeEigenvalueTISE(MathLib& mathlib) : 
    _MathLib(mathlib) {
}
void GeneralizeEigenvalueTISE::Solve() {
    double memory = 3.*(2*_order-1)*_N + 2.;
    memory = memory*16/1024./1024./1024.;
    Log::info("Estimated memory required: " + std::to_string(memory) + " GB.");

    Matrix H0 = _MathLib.CreateMatrix(_N, _N, 2*_order-1);
    Matrix S = _MathLib.CreateMatrix(_N, _N, 2*_order-1);
    Matrix temp = _MathLib.CreateMatrix(_N, _N, 2*_order-1);

    // OPTIMIZATION: these are banded
    std::vector<complex> kinBlockStore(_N*_N);
    std::vector<complex> r2BlockStore(_N*_N);

    // laplacian part is always the same (only depends on i,j) so cache
    for (int i = 0; i < _N; i++) {
        for (int j = i; j < _N; j++) {
            kinBlockStore[i + j*_N] = kinBlockStore[j + i*_N] = _basis.Integrate(i+1, j+1, 1,1) / 2.0;
            r2BlockStore [i + j*_N] =  r2BlockStore[j + i*_N] = _basis.Integrate(i+1, j+1, [] (complex r) {
                return 1./r/r;
            });
        }
    }

    // fill overlap matrix
    Log::info("Building overlap matrix.");
    S->FillBandedBlock(_order-1, _N, [=](int row, int col) {
        int i, j, l1, l2, m1, m2;
        // ILMFrom(row, i, l1, m1);                // used for non-central potential
        // ILMFrom(col, j, l2, m2);                // used for non-central potential
        i = row;
        j = col;
        
        return _basis.Integrate(i+1, j+1);
    }); //, _N, _lmax, _mmax);

    
    // MatView(std::dynamic_pointer_cast<PetscMatrix>(S)->_petsc_mat, 0);
    // exit(0);

    _values.resize(_nmax);
    _vectors.resize(_nmax);
    
    // FIX: This assumes a central potential
    // for each l we get several n's
    for (int l = 0; l < _nmax; l++) {
        Log::info("Building Hamiltonian matrix (l=" + std::to_string(l) + ").");
        // ----------------
        // Fill H0 with kinetic energy
        H0->FillBandedBlock(_order-1, _N, [=](int row, int col) {
            int i, j, l1, l2, m1, m2;
            // ILMFrom(row, i, l1, m1);                // used for non-central potential
            // ILMFrom(col, j, l2, m2);                // used for non-central potential
            i = row;
            j = col;
            
            // kinetic energy and centrifugal term
            return kinBlockStore[i + j*_N] + 0.5*l*(l+1.)*r2BlockStore[i + j*_N];
        }); //, _N, _lmax, _mmax);

    // MatView(std::dynamic_pointer_cast<PetscMatrix>(H0)->_petsc_mat, 0);
    // exit(0);
        // now add all the potential terms`1
        for (auto& p : _potentials) {
            p->FillMatrix(_basis, temp, _N);    // get potential matrix
            _MathLib.AXPY(H0, 1., temp);                             // add potential to H0
        }

        
        _MathLib.Eigen(H0, S, _nmax - l, _tol, _values[l], _vectors[l]);
        // for (int n = l+1; n <= _nmax; n++) {
        //     std::cout << "eigen value: (" << n << ", " << l << ") = " << value[l][n-l-1] << std::endl;
        // }

    }
}

#include <iostream>

void GeneralizeEigenvalueTISE::Store() const {
    // Overwrites any existing file
    // maybe unnecessary?
    std::stringstream ss;
    auto hdf5 = _MathLib.OpenHDF5(_output_filename, 'w');
    hdf5->PushGroup("vectors");
    for (int l = 0; l < _nmax; l++) {                   // l-quantum number
        auto& l_block = _vectors[l];
        for (int n = l+1; n <= _nmax; n++) {            // n-quantum number
            int index = n-l-1;

            ss.str("");                                 // clear string stream
            ss << "(" << n << ", " << l << ")";

            hdf5->WriteVector(ss.str().c_str(), l_block[index]);
            hdf5->WriteAttribute(l_block[index], "energy", std::real(_values[l][index]));
        }
    }    
    hdf5->PopGroup();
}


void GeneralizeEigenvalueTISE::Finish() {
    _values.clear();
    _vectors.clear();
}