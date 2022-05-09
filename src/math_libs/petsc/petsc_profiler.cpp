#include "math_libs/petsc/petsc_lib.h"

void PetscProfiler::Push(const std::string& name) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Profiler::Push(name);
}
void PetscProfiler::Pop(const std::string& name) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Profiler::Pop(name);
}
void PetscProfiler::Print() {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Profiler::Print();
}
bool PetscProfiler::PrintTo(const std::string& filename) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Profiler::PrintTo(filename);
}