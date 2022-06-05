#pragma once

#include "maths/maths.h"
#include "bspline/bspline.h"
#include "tdse/potential.h"

class TISE {
protected:
    Basis::BSpline _basis;

    int _N, _order, _nodes;
    int _lmax, _mmax;
    int _lmin;
    bool _expanding;

    double _xmax, _xmin;

    int _nmax;
    double _tol;
    
    std::vector<Potential::Ptr_t> _potentials;

    std::string _output_filename;
public:
    typedef std::shared_ptr<TISE> Ptr_t;

    TISE();

    void SetupBasis(double xmin, double xmax, 
                    int order, int nodes,
                    double ecs_r0, double ecs_theta,
                    Basis::BSpline::Sequence seq);
    virtual void Solve() = 0;

    void AddPotential(Potential::Ptr_t pot);

    void SetTolerance(double tolerance);
    void SetNMax(int nmax);
    void SetLMax(int lmax);
    void SetLMin(int lmin);
    void SetFilename(const std::string& filename);
    void SetExpanding(bool flag);
    
    virtual void Store() const = 0;
    virtual void Finish() = 0;
    
};
