#include "math_libs/petsc/petsc_lib.h"
#include <iostream>

PetscSolver::PetscSolver(int restart_iter, int max_iter) {
    PetscErrorCode ierr;
    ierr = KSPCreate(PETSC_COMM_WORLD,&_petsc_ksp);PETSCASSERT(ierr);
    ierr = KSPSetType(_petsc_ksp, KSPGMRES);PETSCASSERT(ierr);
    ierr = KSPGMRESSetRestart(_petsc_ksp, restart_iter);PETSCASSERT(ierr);
    // ierr = KSPSetPCSide(ksp,PC_SYMMETRIC);CHKERRQ(ierr);
    ierr = KSPGetPC(_petsc_ksp,&_petsc_pc);PETSCASSERT(ierr);
    ierr = PCSetType(_petsc_pc,PCJACOBI);PETSCASSERT(ierr);
    ierr = KSPSetPC(_petsc_ksp,_petsc_pc);PETSCASSERT(ierr);
    // ierr = KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);CHKERRQ(ierr);
    // ierr = KSPSetTolerances(_petsc_ksp,1e-15, 1e-40, 10000., 10000);PETSCASSERT(ierr);
    ierr = KSPSetTolerances(_petsc_ksp, 1.e-15, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);PETSCASSERT(ierr);
    ierr = KSPSetInitialGuessNonzero(_petsc_ksp,PETSC_TRUE);PETSCASSERT(ierr);
    ierr = KSPSetFromOptions(_petsc_ksp);PETSCASSERT(ierr);
}
PetscSolver::~PetscSolver() {
    PetscErrorCode ierr;

    ierr = KSPDestroy(&_petsc_ksp);PETSCASSERT(ierr);
}

void PetscSolver::SetBlockedPC(int blocks) {
    PetscErrorCode ierr;
    ierr = PCSetType(_petsc_pc,PCBJACOBI);PETSCASSERT(ierr);
    ierr = PCBJacobiSetTotalBlocks(_petsc_pc,blocks,NULL);PETSCASSERT(ierr);
    ierr = KSPSetPC(_petsc_ksp,_petsc_pc);PETSCASSERT(ierr);
}

bool PetscSolver::Solve(const Matrix A, const Vector b, Vector x) {
    PetscErrorCode ierr;

    auto petscA = std::dynamic_pointer_cast<PetscMatrix>(A);
    auto petscb = std::dynamic_pointer_cast<PetscVector>(b);
    auto petscx = std::dynamic_pointer_cast<PetscVector>(x);

    ierr = KSPSetOperators(_petsc_ksp, petscA->_petsc_mat, petscA->_petsc_mat);PETSCASSERT(ierr);
    ierr = KSPSolve(_petsc_ksp, petscb->_petsc_vec, petscx->_petsc_vec);PETSCASSERT(ierr);

    KSPGetConvergedReason(_petsc_ksp, &_reason);
    if (_reason < 0) {
        KSPGetConvergedReasonString(_petsc_ksp, &_strreason);
        std::cout << "Divergence Reason: " << _strreason << std::endl;
        return false;
    }
    return true;
}