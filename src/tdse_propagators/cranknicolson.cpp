
#include "cranknicolson.h"
#include <limits>

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include "utility/index_manip.h"
#include "utility/logger.h"
#include "utility/profiler.h"


using namespace std::complex_literals;

CrankNicolsonTDSE::CrankNicolsonTDSE(MathLib& lib) : TDSE(lib) {
}
void CrankNicolsonTDSE::Initialize() {
    ProfilerPush();

    double memory = 2*(2*_order-1) + 4*_maxBands;
    if (_pol[X]) memory += 8*_order-4;
    if (_pol[Y]) memory += 8*_order-4;
    if (_pol[Z]) memory += 4*_order-2;
    memory += 2; 
    memory *= _dof;
    memory = memory*16/1024./1024./1024.;

    Log::info("Initialize TDSE");
    // ---------------------------------------------------------------
    // initialize field free, overlap, and interaction matrices
    Log::info("Building Hamiltonian and overlap matrix...");
    Log::info("Estimated memory required: " + std::to_string(memory) + " GB.");
    Log::info("Allocating space...");
    Matrix H0 = _MathLib.CreateMatrix(_dof, _dof, 2*_order-1);
    Matrix S = _MathLib.CreateMatrix(_dof, _dof, 2*_order-1);
    
    if (_pol[X])
        _HI[X] = _MathLib.CreateMatrix(_dof, _dof, 8*_order-4);
    if (_pol[Y])
        _HI[Y] = _MathLib.CreateMatrix(_dof, _dof, 8*_order-4);
    if (_pol[Z])
        _HI[Z] = _MathLib.CreateMatrix(_dof, _dof, 4*_order-2);
    
    _U0p = _MathLib.CreateMatrix(_dof, _dof, _maxBands);
    _U0m = _MathLib.CreateMatrix(_dof, _dof, _maxBands);
    _Up = _MathLib.CreateMatrix(_dof, _dof, _maxBands);
    _Um = _MathLib.CreateMatrix(_dof, _dof, _maxBands);

    Log::info("...");
    FillFieldFree(H0);
    Log::info("...");
    FillOverlap(S);
    Log::info("...");

    if (_pol[X])
        FillInteractionX(_HI[X]);
    if (_pol[Y])
        FillInteractionY(_HI[Y]);
    if (_pol[Z])
        FillInteractionZ(_HI[Z]);

    // ---------------------------------------------------------------
    // Initialize the static propagator matrices (U0+/-)
    Log::info("Building propagator matrix...");

    FillU0(_U0p);
    FillU0(_U0m);

    _U0p->AssembleBegin();
    _U0p->AssembleEnd();
    _U0m->AssembleBegin();
    _U0m->AssembleEnd();


    // add overlap (and zero off-diagonal blocks)
    _MathLib.AYPX(_U0p, 0, S);
    _MathLib.AYPX(_U0m, 0, S);
    // add fieldfree hamiltonian
    _MathLib.AXPY(_U0p, 0.5i*_dt, H0);
    _MathLib.AXPY(_U0m, -0.5i*_dt, H0);

    _Up->Duplicate(_U0p);
    _Um->Duplicate(_U0m);
    // ---------------------------------------------------------------
    // Load initial state 
    // TODO: possibly move this to base class. All propagators need an initial state
    Log::info("Loading initial state...");
    Vector temp = _MathLib.CreateVector(_N);        // vector from hdf5 file
    std::vector<Vector> vecs(_dof/_N);             // these will all be concatenated
    for (auto& v : vecs) {
        v = _MathLib.CreateVector(_N); 
        v->Zero();
    }

    
    auto hdf5 = _MathLib.OpenHDF5(_initial_state_filename, 'r');
    hdf5->PushGroup("vectors");

    // compute norm from amplitudes
    double norm = 0;
    for (auto& state : _initial_state)
        norm += state.amplitude*state.amplitude;

    // merge the initial eigenstates into 'vecs'
    std::stringstream name_ss;
    for (auto& state : _initial_state) {
        name_ss.str("");                                         // clear string stream
        name_ss << "(" << state.n << ", " << state.l << ")";     // name of state

        int mBlock = RowFrom(state.m, _Ms, _mRows)/_N;
        int lBlock = mBlock + (state.l-std::abs(state.m));
        hdf5->ReadVector(name_ss.str().c_str(), temp);                  // read in the state
        temp->Scale(state.amplitude*std::exp(1.i*state.phase));         // scale by norm, amplitude and phase
        _MathLib.AXPY(vecs[lBlock], 1., temp);                          // sum with other similar l's
    }
    
    hdf5->PopGroup();

    // concatenate into "_psi"
    _psi = _MathLib.CreateVector(_dof);             // _dof needs to be the sum of all the vectors in 'vecs'
    _psi->Concatenate(vecs);                        // append all the initial vectors together
    _psi->Scale(1./sqrt(norm));                     // and normalize
    _psi_temp = _MathLib.CreateVector(_dof);        // storage used to hold intermediate psi during propagation
    //-----------------------------------------------
    // Create solver
    _solver = _MathLib.CreateGMRESSolver();
    _solver->SetBlockedPC(_N);
    //-----------------------------------------------
    Log::info("Crank-Nicolson initialization complete.");

    // z -> m=m, l=l+-1
    // x -> m=m+-1, l=l+-1
    // y -> m=m+-1, l=l+-1
    
    ProfilerPop();
}

void CrankNicolsonTDSE::Finish() {
    _U0p = nullptr;
    _U0m = nullptr;
    _HI[X] = nullptr;
    _HI[Y] = nullptr;
    _HI[Z] = nullptr;
    _Up = nullptr;
    _Um = nullptr;
    _psi = nullptr;
    _psi_temp = nullptr;
    _solver = nullptr;
}
bool CrankNicolsonTDSE::DoStep(int it, double t, double dt) {
    _Up->Copy(_U0p);
    _Um->Copy(_U0m);

    for (int xn = X; xn <= Z; xn++) {
        if (_HI[xn]) {
            
            _MathLib.AXPY(_Up, 0.5i*dt*(-1.i*_field[xn][it]), _HI[xn]);
            _MathLib.AXPY(_Um, -0.5i*dt*(-1.i*_field[xn][it]), _HI[xn]);
        }
    }
    
    _MathLib.Mult(_Um, _psi, _psi_temp);
    if (!_solver->Solve(_Up, _psi_temp, _psi)) {
        std::cout << "divergence!" << std::endl;
        return false;               // failure
    }



    return true;
}

