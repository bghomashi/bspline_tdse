
#include "eigen_solvers/gen_eigen_tise.h"
#include "utility/logger.h"
#include "utility/profiler.h"
#include "utility/file_exists.h"

#include <string>
#include <sstream>
#include <complex>
#include <iomanip>

GeneralizeEigenvalueTISE::GeneralizeEigenvalueTISE(MathLib& mathlib) : 
    _MathLib(mathlib) {
}
void GeneralizeEigenvalueTISE::Solve() {
    // dump the whole field into a text file
    {
        int numGrid = 500;
        ASCII txt_file = ASCII(_MathLib.OpenASCII("tise_basis.txt", 'w'));
        int num_bsplines = _basis.getNumBSplines();
        std::vector<double> grid(numGrid);
        std::vector<std::vector<complex>> splines(num_bsplines);

        double dx = 30.0/(numGrid - 1);
        for (int i = 0; i < numGrid; i++)
            grid[i] = i*dx;
            
        for (int bs = 0; bs < splines.size(); bs++)
            splines[bs] = _basis.getBSpline(grid, bs+1);

        // set up stringstream for formating and grab some shortcut to values we need
        std::stringstream ss;
        ss << std::setprecision(8) << std::scientific;
        
        // loop though the entire pulse
        for (int i = 0; i < numGrid; i++) {
            ss.str("");
            ss << grid[i];
            for (auto& bs : splines) {
                ss << "\t" << std::real(bs[i]);
                ss << "\t" << std::imag(bs[i]);
            }
            ss << "\n";
            txt_file->Write(ss.str());
        }

        // done
        txt_file = nullptr;
    }

    ProfilerPush();

    double memory = 3.*(2*_order-1)*_N + 2.;
    memory = memory*16/1024./1024./1024.;
    LOG_INFO("Estimated memory required: " + std::to_string(memory) + " GB.");

    Matrix H0 = _MathLib.CreateMatrix(_N, _N, 2*_order-1);
    Matrix S = _MathLib.CreateMatrix(_N, _N, 2*_order-1);
    Matrix temp = _MathLib.CreateMatrix(_N, _N, 2*_order-1);

    std::vector<int> Ls(std::min(_nmax, _lmax) - _lmin + 1);
    for (int l = _lmin; l <= std::min(_nmax, _lmax); l++)
        Ls[l - _lmin] = l;

    if (!file_exists(_output_filename))
        _expanding = false;
 
    // if we *are* going to expand the basis,
    // don't recompute the L-blocks if we don't need to.
    if (_expanding) {
        std::stringstream ss;
        auto hdf5 = _MathLib.OpenHDF5(_output_filename, 'r');

        for (auto it = Ls.begin(); it != Ls.end();) {           // for each of the L's we want
            int l = *it;
            bool skip_this_L = true; 
            // check if we have enough N's
            for (int n = l+1; n <= _nmax; n++) {                // n-quantum number
                ss.str("");                                     // clear string stream
                ss << "(" << n << ", " << l << ")";

                if (!hdf5->HasVector(ss.str())) {               // if we need an (n,l) that is not found
                    skip_this_L = false;                        // do *not* skip this L
                    break;
                }
            }
            if (skip_this_L)
                it = Ls.erase(it);
            else 
                it++;
        }
    }


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
    LOG_INFO("Building overlap matrix.");
    S->FillBandedBlock(_order-1, _N, [=](int row, int col) {
        int i, j, l1, l2, m1, m2;
        // ILMFrom(row, i, l1, m1);                // used for non-central potential
        // ILMFrom(col, j, l2, m2);                // used for non-central potential
        i = row;
        j = col;
        
        return _basis.Integrate(i+1, j+1);
    });

    _values.resize(_nmax);
    _vectors.resize(_nmax);
    

    // open HDF5 file
    std::stringstream ss;
    
    auto hdf5 = _MathLib.OpenHDF5(_output_filename, (_expanding ? 'a' : 'w'));
    hdf5->PushGroup("vectors");


    // FIX: This assumes a central potential
    // for each l we get several n's
    for (int l : Ls) {
        LOG_INFO("Building Hamiltonian matrix (l=" + std::to_string(l) + ").");
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
        });


        // now add all the potential terms`1
        for (auto& p : _potentials) {
            p->FillMatrix(_basis, temp, _N);    // get potential matrix
            _MathLib.AXPY(H0, 1., temp);                             // add potential to H0
        }

        _MathLib.Eigen(H0, S, _nmax - l, _tol, _values[l], _vectors[l]);
        // output eigen values
        for (int j = 0; j < _values[l].size(); j++) {
            ss.str("");  
            ss << "energy " << j << " : (" << std::real(_values[l][j]) << " ," << std::imag(_values[l][j]) << ")";
            LOG_INFO(ss.str());
        }

        // store these vectors
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
    ProfilerPop();
}

void GeneralizeEigenvalueTISE::Store() const {
    ProfilerPush();
    // Overwrites any existing file
    // maybe unnecessary?
    for (int l = 0; l < _nmax; l++) {                   // l-quantum number
    }    
    
    ProfilerPop();
}


void GeneralizeEigenvalueTISE::Finish() {
    _values.clear();
    _vectors.clear();
}