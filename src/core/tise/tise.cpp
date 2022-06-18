#include "tise/tise.h"

#include <iostream>

TISE::TISE() : _expanding(false), _nmax(1), _tol(1e-10), _lmax(0), _lmin(0) {}
void TISE::SetupBasis(double xmin, double xmax, 
                    int order, int nodes,
                    double ecs_r0, double ecs_theta,
                    Basis::BSpline::Sequence seq,
                    double seq_parameter) {
    _xmin = xmin; _xmax = xmax;

    _order = order;
    _nodes = nodes;

    _basis.Initialize(_order, _nodes, 
                      _xmin, _xmax, 
                      seq, 
                      Basis::BSpline::ECS{ecs_r0, ecs_theta},
                      seq_parameter);                           // only used if (seq != linear)
    _basis.setSkipFirst();
	_basis.setSkipLast();
    _N = _basis.getNumBSplines();
}


void TISE::AddPotential(Potential::Ptr_t pot) {
    _potentials.push_back(pot);
}
void TISE::SetTolerance(double tolerance) {
    _tol = tolerance;
}
void TISE::SetNMax(int nmax) {
    _nmax = nmax;
}
void TISE::SetLMax(int lmax) {
    _lmax = lmax;
}
void TISE::SetLMin(int lmin) {
    _lmin = lmin;
}
void TISE::SetFilename(const std::string& filename) {
    _output_filename = filename;
}
void TISE::SetExpanding(bool flag) {
    _expanding = flag;
}