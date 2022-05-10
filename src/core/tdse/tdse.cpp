#include "tdse.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "utility/index_manip.h"
#include "utility/logger.h"
#include "utility/profiler.h"

TDSE::TDSE(MathLib& lib) : _MathLib(lib), _cylindricalSymmetry(true), _checkpoints(0) {
    _pol[X] = _pol[Y] = _pol[Z] = false;
}
void TDSE::SetupBasis(double xmin, double xmax, 
                    int order, int nodes, 
                    int lmax, int mmax, 
                    double ecs_r0, double ecs_theta,
                    Basis::BSpline::Sequence seq) {
    _xmin = xmin; _xmax = xmax;

    _order = order;
    _nodes = nodes;
    _lmax = lmax;
    _mmax = mmax;

    _basis.Initialize(_order, _nodes, 
                      _xmin, _xmax, 
                      seq, 
                      Basis::BSpline::ECS{ecs_r0, ecs_theta});
    _basis.setSkipFirst();
    _basis.setSkipLast();
    _N = _basis.getNumBSplines();

    // first check the cylindrical symmetry is broken.
    // - we *could* rotate any one axis to the z-axis to preserve symmetry
    // - this assumes the potential is at least cylindrically symmetric also
    for (const auto& p : _pulses) {
        if (p->polarization_vector[X]) _pol[X] = true;
        if (p->polarization_vector[Y]) _pol[Y] = true;
        if (p->polarization_vector[Z]) _pol[Z] = true;
    }
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
const std::vector<int>& TDSE::Ms() const {
    return _Ms;
}
const std::vector<int>& TDSE::MRows() const {
    return _mRows;
}
const bool* TDSE::Polarization() const {
    return _pol;
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
const std::string& TDSE::GetInitialStateFile() const {
    return _initial_state_filename;
}
int TDSE::GetInitialStateNmax() const {
    return _eigen_state_nmax;
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
    Log::info("Beginning propagation...");

    double t = 0;
    int NT;
    // find total length of simulation by finding the latest nonzero pulse.
    _tmin = _tmax = 0; 
    for (auto& p : _pulses)
        _tmax = std::max(_tmax, p->delay + p->duration);
    NT =  (_tmax - _tmin)/ _dt + 1;

    Log::info("Calculating cummulative field...");
    // calculate the field for each timestep in advance
    for (auto& p : _pulses) {
        if (p->polarization_vector[X] != 0.)
            _field[X].resize(NT);
        if (p->polarization_vector[Y] != 0.)
            _field[Y].resize(NT);
        if (p->polarization_vector[Z] != 0.)
            _field[Z].resize(NT);
    }

    //std::ofstream file("data/x_field.txt");
    //file << std::setprecision(8) << std::scientific;
    if (_field[X].size() > 0) {
        for (int it = 0; it < NT; it++) {
            for (auto& p : _pulses)
                _field[X][it] += (*p)(it*_dt)*p->polarization_vector[X];
            //file << _field[X][it] << std::endl;
        }
    }
    if (_field[Y].size() > 0) {
        for (int it = 0; it < NT; it++)
            for (auto& p : _pulses)
                _field[Y][it] += (*p)(it*_dt)*p->polarization_vector[Y];
    }

    if (_field[Z].size() > 0) {
        for (int it = 0; it < NT; it++)
            for (auto& p : _pulses)
                _field[Z][it] += (*p)(it*_dt)*p->polarization_vector[Z];

    }
    //exit(0);

    for (auto& obs : _observables)
        obs->Startup();

    Log::info("Propagation for " + std::to_string(NT) + " timesteps...");
    Profile::Push("Total time stepping");
    
    // do simulation
    for (int it = 0; it < NT; it++) {
        t = it*_dt;
        if (!DoStep(it, t, _dt)) break;
        DoCheckpoint(it, NT);
        DoObservables(it, t, _dt);
    }
    
    Profile::Pop("Total time stepping");

    for (auto& obs : _observables)
        obs->Shutdown();

    Finish();
}
void TDSE::DoCheckpoint(int it, int NT) {
    if ((_checkpoints != 0) && (it % _checkpoints == 0)) {
        Log::info("iteration: " + std::to_string(it) + "/" + std::to_string(NT));
        for (auto& obs : _observables)
            obs->Flush();
    }
}
void TDSE::DoObservables(int it, double t, double dt) {
    for (auto& obs : _observables)
        obs->DoObservable(it, t, dt);
}
