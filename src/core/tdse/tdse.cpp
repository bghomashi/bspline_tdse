#include "tdse.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <complex>
#include "utility/index_manip.h"
#include "utility/logger.h"
#include "utility/profiler.h"
#include "utility/file_exists.h"

#include "math_libs/petsc/petsc_lib.h"

using namespace std::complex_literals;


TDSE::TDSE(MathLib& lib) : _MathLib(lib), _do_propagate(true), _restarting(false), _cylindricalSymmetry(true), _checkpoints(0) {
    _pol[X] = _pol[Y] = _pol[Z] = false;
    _ecs_r0 = 0;
    _ecs_theta = 0;
}
void TDSE::SetupBasis(double xmin, double xmax, 
                    int order, int nodes, 
                    int lmax, int mmax, 
                    double ecs_r0, double ecs_theta,
                    Basis::BSpline::Sequence seq,
                    double seq_parameter) {
    _xmin = xmin; _xmax = xmax;

    _order = order;
    _nodes = nodes;
    _lmax = lmax;
    _mmax = mmax;

    _basis.Initialize(_order, _nodes, 
                      _xmin, _xmax, 
                      seq, 
                      Basis::BSpline::ECS{ecs_r0, ecs_theta},
                      seq_parameter);
    _basis.setSkipFirst();
    _basis.setSkipLast();
    _N = _basis.getNumBSplines();

    // first check the cylindrical symmetry is broken.
    // - we *could* rotate any one axis to the z-axis to preserve symmetry
    // - this assumes the potential is at least cylindrically symmetric also
    for (const auto& p : _pulses) {
        if (p->polarization_vector.x || p->minor_polarization_vector.x) _pol[X] = true;
        if (p->polarization_vector.y || p->minor_polarization_vector.y) _pol[Y] = true;
        if (p->polarization_vector.z || p->minor_polarization_vector.z) _pol[Z] = true;
    }
    // if (_pol[X] || _pol[Y]) {
    //     LOG_CRITICAL("DOING 3D calculation");
    //     exit(-1);
    // }
    // if there is any polarization in the x/y direction
    if (_pol[X] || _pol[Y])
        _cylindricalSymmetry = false;

    // if the interaction preserves cylindrical symmetry
    if (_cylindricalSymmetry) {
        // the Ms only come from the initial state and are not coupled
        for (auto& state : _initial_state)
            _Ms.push_back(state.m);
    
        // sort and remove duplicates
        std::sort(_Ms.begin(), _Ms.end());
        _Ms.erase( std::unique(_Ms.begin(), _Ms.end() ), _Ms.end() );

        // still each m-block is independent so each row of propagator
        // has 6*order - 3 elements
        _maxBands = 6*_order - 3;
    } else {                // no symmetry - either due to potential or laser
        // we need all the m-values anyway
        _Ms.reserve(2*mmax+1);
        for (int m = -mmax; m <= mmax; m++)
            _Ms.push_back(m);

        
        // every m-block is coupled the one before and after it
        _maxBands = 12*_order - 6;
    }

    // table to quickly find the first row for an m-block
    _mRows.reserve(_Ms.size());

    // count degrees of freedom and 
    // generate row table for m's
    _dof = 0;
    for (const int m : _Ms) {
        _mRows.push_back(_dof);
        _dof += _N*(_lmax - std::abs(m) + 1);
    }
}

const Vector TDSE::Psi() const {
    return _psi;
}
MathLib& TDSE::MathLibrary() {
    return _MathLib;
}
const Basis::BSpline& TDSE::Basis() const {
    return _basis;
}
int TDSE::DOF() const {
    return _dof;
}
int TDSE::Lmax() const {
    return _lmax;
}
int TDSE::Mmax() const {
    return _mmax;
} 
double TDSE::Xmin() const {
    return _xmin;
}
double TDSE::Xmax() const {
    return _xmax;
}
double TDSE::Tmin() const {
    return _tmin;
}
double TDSE::Tmax() const {
    return _tmax;
}
double TDSE::Timestep() const {
    return _dt;
}
const std::vector<int>& TDSE::Ms() const {
    return _Ms;
}
const std::vector<int>& TDSE::MRows() const {
    return _mRows;
}
const bool* TDSE::Polarization() const {
    return _pol;
}
const std::vector<double>& TDSE::GetField(int dim_index) const {
    return _field[dim_index];
}

int TDSE::NumTimeSteps() const {
    return _NT;
}
const std::vector<Pulse::Ptr_t>& TDSE::Pulses() const {
    return _pulses;
}

void TDSE::AddPulse(Pulse::Ptr_t p) {
    _pulses.push_back(p);
}
void TDSE::AddInitialState(int n, int l, int m, double phase,  double amplitude) {
    _initial_state.emplace_back(state_descriptor{n, l, m, phase, amplitude});
}
void TDSE::AddObservable(Observable::Ptr_t obs) {
    _observables.push_back(obs);
}
void TDSE::AddPotential(Potential::Ptr_t pot) {
    _potentials.push_back(pot);
}
void TDSE::SetInitialStateFile(const std::string& filename) {
    _initial_state_filename = filename;
}
void TDSE::SetEigenStateNmax(int nmax) {
    _eigen_state_nmax = nmax;
}
void TDSE::SetEigenStateLmax(int lmax) {
    _eigen_state_lmax = lmax;
}
const std::string& TDSE::GetInitialStateFile() const {
    return _initial_state_filename;
}
int TDSE::GetInitialStateNmax() const {
    return _eigen_state_nmax;
}
int TDSE::GetInitialStateLmax() const {
    return _eigen_state_lmax;
}
void TDSE::SetRestart(bool flag) {
    _restarting = flag;
}
void TDSE::SetDoPropagate(bool flag) {
    _do_propagate = flag;
}
void TDSE::SetECS(double ecs_r0, double ecs_theta) {
    _ecs_r0 = ecs_r0;
    _ecs_theta = ecs_theta;
}
void TDSE::SetTimestep(double dt) {
    _dt = dt;
}
void TDSE::SetCheckpoints(int checkpoint) {
    _checkpoints = checkpoint;
}
const std::vector<Potential::Ptr_t>& TDSE::Potentials() const {
    return _potentials;
}

void TDSE::Propagate() {
    double t = 0;
    int start_iteration = 0;

    // find total length of simulation by finding the latest nonzero pulse.
    _tmin = _tmax = 0; 
    for (auto& p : _pulses)
        _tmax = std::max(_tmax, p->delay + p->duration);
    _NT =  (_tmax - _tmin)/ _dt + 1;

    _psi = _MathLib.CreateVector(_dof);

    // are we restarting?
    if (_restarting)
        _restarting = file_exists("TDSE.h5");
    
    _tdse_out = _MathLib.OpenHDF5("TDSE.h5", (_restarting ? 'a' : 'w'));

    if (_restarting) {
        if (!CompareTDSEH5wInput())
            exit(-1);           // failure
    } else 
        WriteParametersToTDSE();


    if (!_restarting || !LoadLastCheckpoint(start_iteration)) {
        LoadInitialState();
        WriteInitialState();
    }


    LOG_INFO("Beginning propagation...");


    // calculate the field for each timestep in advance
    LOG_INFO("Calculating cummulative field...");
    ComputeFields();

    // allow observables to initialize
    for (auto& obs : _observables)
        obs->Startup(start_iteration);

    if (!_do_propagate) {
        _tdse_out = nullptr;

        // allow observables to complete
        for (auto& obs : _observables)
            obs->Shutdown();
        _tdse_out = nullptr;

        Finish();
        return;
    }

    // begin
    LOG_INFO("Propagation for " + std::to_string(_NT) + " timesteps...");
    Profile::Push("Total time stepping");
    
    // do simulation
    for (int it = start_iteration; it < _NT; it++) {
        t = it*_dt;
        if (!DoStep(it, t, _dt)) break;
        DoCheckpoint(it);
        DoObservables(it, t, _dt);
    }
    
    Profile::Pop("Total time stepping");

    WriteFinalState();

    _tdse_out = nullptr;

    // allow observables to complete
    for (auto& obs : _observables)
        obs->Shutdown();

    // finish up
    Finish();
}
void TDSE::DoCheckpoint(int it) {
    if ((_checkpoints != 0) && (it % _checkpoints == 0)) {
        LOG_INFO("iteration: " + std::to_string(it) + "/" + std::to_string(_NT));
        for (auto& obs : _observables)
            obs->Flush();                           // dumb all the observables

        // dump psi to file
        _tdse_out->PushGroup("checkpoints");
        _tdse_out->WriteVector(std::to_string(it), _psi);
        _tdse_out->PopGroup();

        // write last checkpoint iteration
        _tdse_out->PushGroup("parameters");
        _tdse_out->WriteAttribute("last_checkpoint", it);           // should update
        _tdse_out->PopGroup();
    }
}
void TDSE::DoObservables(int it, double t, double dt) {
    for (auto& obs : _observables)
        obs->DoObservable(it, t, dt);
}

void TDSE::ComputeFields() {
    for (auto& p : _pulses) {
        if (p->polarization_vector.x != 0. || p->minor_polarization_vector.x != 0.)
            _field[X].resize(_NT);
        if (p->polarization_vector.y != 0. || p->minor_polarization_vector.y != 0.)
            _field[Y].resize(_NT);
        if (p->polarization_vector.z != 0. || p->minor_polarization_vector.z != 0.)
            _field[Z].resize(_NT);
    }

    // if (_field[X].size() > 0 || _field[Y].size() > 0) {
    //     LOG_CRITICAL("DOING 3D calculation");
    //     exit(-1);
    // }

    if (_field[X].size() > 0) {
        for (int it = 0; it < _NT; it++) {
            for (auto& p : _pulses)
                _field[X][it] += (*p)(it*_dt).x;
        }
    }
    if (_field[Y].size() > 0) {
        for (int it = 0; it < _NT; it++)
            for (auto& p : _pulses)
                _field[Y][it] += (*p)(it*_dt).y;
    }
    if (_field[Z].size() > 0) {
        for (int it = 0; it < _NT; it++)
            for (auto& p : _pulses)
                _field[Z][it] += (*p)(it*_dt).z;
    }
}



void TDSE::LoadInitialState() {
    LOG_INFO("Loading initial state...");
    if (!file_exists(_initial_state_filename)) {
        LOG_CRITICAL("eigenstate file does not exists: " + _initial_state_filename);
        exit(-1);
    }

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
        temp->Scale(state.amplitude*std::exp(1.i*state.phase));         // scale by amplitude and phase
        _MathLib.AXPY(vecs[lBlock], 1., temp);                          // sum with other similar l's
    }
    hdf5->PopGroup();

    // concatenate into "_psi"
    _psi->Concatenate(vecs);                        // append all the initial vectors together
    _psi->Scale(1./sqrt(norm));                     // and normalize
}
bool TDSE::LoadLastCheckpoint(int& start_iteration) {
    int it;
    LOG_INFO("Loading last checkpoint...");
    // what was the last successful checkpoint?
    _tdse_out->PushGroup("parameters");
    _tdse_out->ReadAttribute("last_checkpoint", &it);
    _tdse_out->PopGroup();

    if (it == -1) {             // ran but never had a checkpoint
        LOG_INFO("No checkpoints found...");
        return false;
    }
    start_iteration = it;
    // load that checkpoint
    _tdse_out->PushGroup("checkpoints");
    _tdse_out->ReadVector(std::to_string(start_iteration), _psi);
    _tdse_out->PopGroup();

    return true;
}

void TDSE::WriteInitialState() const {
    _tdse_out->PushGroup("initial_state");
    _tdse_out->WriteVector("wavefunction", _psi);
    _tdse_out->PopGroup();
}
void TDSE::WriteFinalState() const {
    _tdse_out->PushGroup("final_state");
    _tdse_out->WriteVector("wavefunction", _psi);
    _tdse_out->PopGroup();
}



void TDSE::WriteParametersToTDSE() const {
    _tdse_out->PushGroup("parameters");
    _tdse_out->WriteAttribute("time_step", _dt);
    _tdse_out->WriteAttribute("time_min", _tmin);
    _tdse_out->WriteAttribute("time_max", _tmax);
    _tdse_out->WriteAttribute("num_timesteps", _NT);
    _tdse_out->WriteAttribute("x_min", _xmin);
    _tdse_out->WriteAttribute("x_max", _xmax);
    _tdse_out->WriteAttribute("nodes", _nodes);
    _tdse_out->WriteAttribute("order", _order);
    _tdse_out->WriteAttribute("l_max", _lmax);
    _tdse_out->WriteAttribute("m_max", _mmax);
    _tdse_out->WriteAttribute("ecs_r0", _ecs_r0);
    _tdse_out->WriteAttribute("ecs_theta", _ecs_theta);
    _tdse_out->WriteAttribute("last_checkpoint", -1);        // for later
    _tdse_out->PopGroup();
}


#define COMPARE_PARAM(x,y) if (x != y) { \
        LOG_CRITICAL("input "#x" does not match TDSE.h5. Expected " + std::to_string(x)); \
        return false; \
    }

bool TDSE::CompareTDSEH5wInput() const {
    // read parameters and compare to input file
    double ecs_r0, ecs_theta;
    int order, nodes;
    int lmax, mmax;
    double xmax, xmin;
    double dt, tmin, tmax;
    int NT;

    _tdse_out->PushGroup("parameters");
    _tdse_out->ReadAttribute("time_step", &dt);
    _tdse_out->ReadAttribute("time_min", &tmin);
    _tdse_out->ReadAttribute("time_max", &tmax);
    _tdse_out->ReadAttribute("num_timesteps", &NT);
    _tdse_out->ReadAttribute("x_min", &xmin);
    _tdse_out->ReadAttribute("x_max", &xmax);
    _tdse_out->ReadAttribute("nodes", &nodes);
    _tdse_out->ReadAttribute("order", &order);
    _tdse_out->ReadAttribute("l_max", &lmax);
    _tdse_out->ReadAttribute("m_max", &mmax);
    _tdse_out->ReadAttribute("ecs_r0", &ecs_r0);
    _tdse_out->ReadAttribute("ecs_theta", &ecs_theta);
    _tdse_out->PopGroup();

    // make sure everything matches
    COMPARE_PARAM(dt,_dt);
    COMPARE_PARAM(tmin,_tmin);
    COMPARE_PARAM(tmax,_tmax);
    COMPARE_PARAM(NT,_NT);
    COMPARE_PARAM(xmin,_xmin);
    COMPARE_PARAM(xmax,_xmax);
    COMPARE_PARAM(nodes,_nodes);
    COMPARE_PARAM(order,_order);
    COMPARE_PARAM(lmax,_lmax);
    COMPARE_PARAM(mmax,_mmax);
    COMPARE_PARAM(ecs_r0,_ecs_r0);
    COMPARE_PARAM(ecs_theta,_ecs_theta);
    
    return true;
}
