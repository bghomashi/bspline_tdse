#include "math_libs/petsc/petsc_lib.h"


PetscVector::PetscVector() {
    _len = 0;
    _petsc_vec = 0;
    _petsc_ctx = 0;
    _petsc_sca_vec =0;
    _petsc_is = 0;
}
PetscVector::PetscVector(int length) {
    PetscErrorCode ierr;
    _petsc_ctx = 0;
    _petsc_sca_vec =0;
    ierr = VecCreate(PETSC_COMM_WORLD,&_petsc_vec);PETSCASSERT(ierr);
    ierr = VecSetSizes(_petsc_vec,PETSC_DECIDE,length);PETSCASSERT(ierr);
    ierr = VecSetFromOptions(_petsc_vec);PETSCASSERT(ierr);
    ierr = VecSet(_petsc_vec, 0.0);PETSCASSERT(ierr);
    ierr = VecAssemblyBegin(_petsc_vec); PETSCASSERT(ierr);
    ierr = VecAssemblyEnd(_petsc_vec); PETSCASSERT(ierr);
    _len = length;
}
PetscVector::~PetscVector() {
    PetscErrorCode ierr;
    ierr = VecScatterDestroy(&_petsc_ctx); PETSCASSERT(ierr);
    ierr = VecDestroy(&_petsc_sca_vec); PETSCASSERT(ierr);
    ierr = VecDestroy(&_petsc_vec); PETSCASSERT(ierr);
}

complex PetscVector::Get(int index) const {
    return 0.;
}
void PetscVector::Get(std::vector<complex>& out) {
    if (_petsc_sca_vec == 0)

    out.resize(_len);

}
void PetscVector::Set(int index, complex value) {
    PetscErrorCode ierr = VecSetValue(_petsc_vec, index,  value, INSERT_VALUES); PETSCASSERT(ierr);
}
void PetscVector::Scale(PetscScalar a) {
    VecScale(_petsc_vec, a);
}
void PetscVector::Duplicate(const Vector& o) {
    // destroy old vector if there was one
    PetscErrorCode ierr;
    ierr = VecScatterDestroy(&_petsc_ctx); PETSCASSERT(ierr);
    ierr = VecDestroy(&_petsc_vec); PETSCASSERT(ierr);
    _petsc_ctx = 0; _petsc_vec = 0;

    auto from = std::dynamic_pointer_cast<PetscVector>(o);
    ierr = VecDuplicate(from->_petsc_vec, &_petsc_vec); PETSCASSERT(ierr);
    _len = from->_len;
}
void PetscVector::Copy(const Vector& o) {
    auto from = std::dynamic_pointer_cast<PetscVector>(o);
    PetscErrorCode ierr = VecCopy(from->_petsc_vec,_petsc_vec); PETSCASSERT(ierr);
}
void PetscVector::Zero() {
    PetscErrorCode ierr = VecZeroEntries(_petsc_vec); PETSCASSERT(ierr);
}
void PetscVector::Concatenate(const std::vector<Vector>& vecs) {
    PetscErrorCode ierr;
    int nx = vecs.size();
    // destroy old vector...
    ierr = VecScatterDestroy(&_petsc_ctx); PETSCASSERT(ierr);
    ierr = VecDestroy(&_petsc_vec); PETSCASSERT(ierr);
    _len = 0; _petsc_ctx = 0; _petsc_vec = 0;


    std::vector<Vec> petsc_vec_array(nx);
    for (int i = 0; i < nx; i++) {
        auto v = std::dynamic_pointer_cast<PetscVector>(vecs[i]);
        petsc_vec_array[i] = v->_petsc_vec;
        _len += v->_len;
    }
    ierr = VecConcatenate(nx, petsc_vec_array.data(), &_petsc_vec, NULL); PETSCASSERT(ierr);
}
void PetscVector::CopyTo(std::vector<complex>& values) {
    PetscMPIInt rank;
    values.resize(_len);

    if (_petsc_ctx == 0) 
        CreateScatter();

    Scatter();

    
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); PETSCASSERT(ierr);
    if (rank == 0) {
        complex* ptr;
        ScatterGetArray(&ptr);

        for (int i = 0; i < _len; i++)
            values[i] = ptr[i];

        ScatterRestoreArray(&ptr);
    }
}

#include<iostream>
void PetscVector::Transform(Vector& out, std::function<std::vector<complex>(const std::vector<complex>&)> f) {
    PetscMPIInt rank;
    std::vector<complex> in_values(Length());
    std::vector<complex> out_values(out->Length());

    if (_petsc_ctx == 0) 
        CreateScatter();

    Scatter();

    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); PETSCASSERT(ierr);
    if (rank == 0) {
        complex* ptr;
        ScatterGetArray(&ptr);

        for (int i = 0; i < _len; i++)
            in_values[i] = ptr[i];

        ScatterRestoreArray(&ptr);

        out_values = f(in_values);

        for (int i = 0; i < out_values.size(); i++)
            out->Set(i, out_values[i]);

        out->AssembleBegin();
        out->AssembleEnd();
    }
}

void PetscVector::CreateScatter() {
    VecScatterCreateToZero(_petsc_vec, &_petsc_ctx, &_petsc_sca_vec);
}
void PetscVector::Scatter() {
    VecScatterBegin(_petsc_ctx, _petsc_vec, _petsc_sca_vec, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(_petsc_ctx, _petsc_vec, _petsc_sca_vec, INSERT_VALUES, SCATTER_FORWARD);
}
void PetscVector::ScatterGetArray(PetscScalar** ptr) {
    VecGetArray(_petsc_sca_vec, ptr);
}
void PetscVector::ScatterRestoreArray(PetscScalar** ptr) {
    VecRestoreArray(_petsc_sca_vec, ptr);
}

Vector PetscVector::GetSubVector(int start, int end) {
    PetscVector* result = new PetscVector();
    PetscInt low, high;

    result->_len = end - start;

    VecCreate(PETSC_COMM_WORLD,&result->_petsc_vec);
    VecSetSizes(result->_petsc_vec,PETSC_DECIDE,end - start);	
    VecSetUp(result->_petsc_vec);
    VecGetOwnershipRange(result->_petsc_vec, &low, &high);
    VecDestroy(&result->_petsc_vec);
    

    ISCreateStride(PETSC_COMM_WORLD, high-low, start+low, 1, &result->_petsc_is);
    VecGetSubVector(_petsc_vec, result->_petsc_is, &result->_petsc_vec);

    return Vector(result);
}
void PetscVector::RestoreSubVector(Vector sub) {
    auto petsc_sub = std::dynamic_pointer_cast<PetscVector>(sub);
    VecRestoreSubVector(_petsc_vec, petsc_sub->_petsc_is, &petsc_sub->_petsc_vec);
    ISDestroy(&petsc_sub->_petsc_is);
    petsc_sub->_petsc_vec = 0;
    petsc_sub->_petsc_is = 0;
}


void PetscVector::AssembleBegin() {
    PetscErrorCode ierr;
    ierr = VecAssemblyBegin(_petsc_vec); PETSCASSERT(ierr);
}
void PetscVector::AssembleEnd() {
    PetscErrorCode ierr;
    ierr = VecAssemblyEnd(_petsc_vec); PETSCASSERT(ierr);
}
