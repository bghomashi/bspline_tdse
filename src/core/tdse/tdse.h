#pragma once

#include "bspline/bspline.h"
#include "maths/maths.h"
#include "potential.h"
#include "laser.h"
#include "observable.h"

class TDSE {
protected:
    // an object to hold the initial state info
    struct state_descriptor {
        int n, l, m;
        double phase, amplitude;
    };

    // math library
    MathLib& _MathLib;
    // basis 
    Basis::BSpline _basis;
    bool _cylindricalSymmetry;
    // dimension information
    int _N, _order, _nodes;
    int _lmax, _mmax;
    int _dof, _maxBands;
    std::vector<int> _Ms;
    std::vector<int> _mRows;                    // starting row for each m


    // the grid domain
    double _xmax, _xmin;

    // the time domain
    double _dt, _tmin, _tmax;

    // physical quantities
    bool _pol[DimIndex::NUM];                   // quick access if there is/is not polarization in x,y,z
    std::vector<double> _field[DimIndex::NUM];
    std::vector<Pulse::Ptr_t> _pulses;
    std::vector<Potential::Ptr_t> _potentials;
    std::vector<state_descriptor> _initial_state;
    Vector _psi;

    // where to load initial state from
    std::string _initial_state_filename;
    int _eigen_state_nmax;

    // a list of observables specified in the input file.
    std::vector<Observable::Ptr_t> _observables;
    int _checkpoints;
public:
    typedef std::shared_ptr<TDSE> Ptr_t;

    TDSE(MathLib& lib);
    void SetupBasis(double xmin, double xmax, 
                    int order, int nodes, 
                    int lmax, int mmax, 
                    double ecs_r0, double ecs_theta,
                    Basis::BSpline::Sequence seq);
    void Propagate();
    void AddPulse(Pulse::Ptr_t p);
    void AddPotential(Potential::Ptr_t pot);
    void AddObservable(Observable::Ptr_t obs);
    void SetTimestep(double dt);
    void SetCheckpoints(int checkpoint);
    void SetInitialStateFile(const std::string& filename);
    void SetEigenStateNmax(int nmax);
    const std::string& GetInitialStateFile() const;
    int GetInitialStateNmax() const;

    void AddInitialState(int n, int l, int m, double phase, double amplitude);

    const Vector Psi() const;
    MathLib& MathLibrary();
    const Basis::BSpline& Basis() const;
    const std::vector<Potential::Ptr_t>& Potentials() const;
    
    int DOF() const;
    int Lmax() const;    
    int Mmax() const;  
    const std::vector<int>& Ms() const;
    const std::vector<int>& MRows() const;  
    double Xmin() const; 
    double Xmax() const;
    const bool* Polarization() const;


    virtual void Initialize() = 0;
    virtual bool DoStep(int it, double t, double dt) = 0;
    void DoCheckpoint(int it, int NT);
    void DoObservables(int it, double t, double dt);
    virtual void Finish() = 0;
};