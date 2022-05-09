#include "math_libs/petsc/petsc_lib.h"


void PetscLogger::info(const std::string& text) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Logger::info(text);
}
void PetscLogger::warn(const std::string& text) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Logger::warn(text);
}
void PetscLogger::critical(const std::string& text) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Logger::critical(text);
}
void PetscLogger::debug(const std::string& text) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Logger::debug(text);
}
void PetscLogger::set_logger_file(const std::string& log_file) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Logger::set_logger_file(log_file);
}