#include "math_libs/petsc/petsc_lib.h"
#include "utility/logger.h"

bool Petsc::Startup(int argc, char **args) {
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc,&args,NULL,"Func Test\n"); 
    Log::info("Initializing PETsc.");
    if (ierr) {
        Log::critical("Failed!");
        return false;
    }
    ierr = SlepcInitialize(&argc,&args,NULL,NULL);if (ierr)
    if (ierr) {
        Log::critical("Failed!");
        return false;
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&_size);CHKERRQ(ierr);

    return true;
}
void Petsc::Shutdown() {
    PetscErrorCode ierr;
    ierr = SlepcFinalize();
    ierr = PetscFinalize();
}



Vector Petsc::CreateVector(int N) {
    return Vector(new PetscVector(N));
}
void Petsc::DestroyVector(Vector& v) {
    v = nullptr;                // If there are other references to v, the object is not destroyed
}

Matrix Petsc::CreateMatrix(int rows, int cols, int numBands) {
    return Matrix(new PetscMatrix(rows, cols, numBands));
}
void Petsc::DestroyMatrix(Matrix& m) {
    m = nullptr;                // If there are other references to m, the object is not destroyed
}
GMRESSolver Petsc::CreateGMRESSolver(int restart_iter, int max_iter) {
    return GMRESSolver(new PetscSolver(restart_iter, max_iter));
}
void Petsc::DestroyGMRESSolver(GMRESSolver& m) {
    m = nullptr;
}

HDF5 Petsc::OpenHDF5(const std::string& filename, char mode) {
    return HDF5(new PetscHDF5(filename, mode));
}
void Petsc::CloseHDF5(HDF5& file) {
    file = nullptr;
}


ASCII Petsc::OpenASCII(const std::string& filename, char mode) {
    return ASCII(new PetscASCII(filename, mode));
}
void Petsc::CloseASCII(ASCII& file) {
    file = nullptr;
}

void Petsc::ParallelPrintf(const std::string& text) {
    if (_rank == 0)
        PetscPrintf(PETSC_COMM_SELF, text.c_str());
}


void Petsc::Mult(const Matrix M, const Vector in, Vector out) {
    auto m = std::dynamic_pointer_cast<PetscMatrix>(M);
    auto a = std::dynamic_pointer_cast<PetscVector>(in);
    auto b = std::dynamic_pointer_cast<PetscVector>(out);
    MatMult(m->_petsc_mat, a->_petsc_vec, b->_petsc_vec);
}
void Petsc::Dot(const Vector a, const Vector b, complex& value) {
    auto aa = std::dynamic_pointer_cast<PetscVector>(a);
    auto bb = std::dynamic_pointer_cast<PetscVector>(b);

    VecDot(aa->_petsc_vec, bb->_petsc_vec, &value);
}
void Petsc::AYPX(Matrix Y, complex a, const Matrix X) {
    auto y = std::dynamic_pointer_cast<PetscMatrix>(Y);
    auto x = std::dynamic_pointer_cast<PetscMatrix>(X);

    MatAYPX(y->_petsc_mat, a, x->_petsc_mat, SUBSET_NONZERO_PATTERN);
}
void Petsc::AXPY(Matrix Y, complex a, const Matrix X) {
    auto y = std::dynamic_pointer_cast<PetscMatrix>(Y);
    auto x = std::dynamic_pointer_cast<PetscMatrix>(X);

    MatAXPY(y->_petsc_mat, a, x->_petsc_mat, SUBSET_NONZERO_PATTERN);
}
void Petsc::AYPX(Vector Y, complex a, const Vector X) {
    auto y = std::dynamic_pointer_cast<PetscVector>(Y);
    auto x = std::dynamic_pointer_cast<PetscVector>(X);

    VecAYPX(y->_petsc_vec, a, x->_petsc_vec);
}
void Petsc::AXPY(Vector Y, complex a, const Vector X) {
    auto y = std::dynamic_pointer_cast<PetscVector>(Y);
    auto x = std::dynamic_pointer_cast<PetscVector>(X);

    VecAXPY(y->_petsc_vec, a, x->_petsc_vec);
}
