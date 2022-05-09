#include "tise/tise.h"

#include <iostream>

TISE::TISE() : _nmax(1), _tol(1e-10) {}
void TISE::SetupBasis(double xmin, double xmax, 
                    int order, int nodes,
                    double ecs_r0, double ecs_theta,
                    Basis::BSpline::Sequence seq) {
    _xmin = xmin; _xmax = xmax;

    _order = order;
    _nodes = nodes;

    _basis.Initialize(_order, _nodes, 
                      _xmin, _xmax, 
                      seq, 
                      Basis::BSpline::ECS{ecs_r0, ecs_theta});
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
void TISE::SetFilename(const std::string& filename) {
    _output_filename = filename;
}