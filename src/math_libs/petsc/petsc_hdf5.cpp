#include "math_libs/petsc/petsc_lib.h"

#include <petscviewerhdf5.h>

PetscHDF5::PetscHDF5() : _viewer(0) {}
PetscHDF5::PetscHDF5(const std::string& filename, char mode) : _viewer(0) {
    PetscErrorCode ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(),
        (mode == 'a' ? FILE_MODE_APPEND :
        (mode == 'w' ? FILE_MODE_WRITE :
                       FILE_MODE_READ))
        , &_viewer); PETSCASSERT(ierr);
    PetscViewerHDF5SetSPOutput(_viewer, PETSC_FALSE);
    PetscViewerSetFromOptions(_viewer);
}
PetscHDF5::~PetscHDF5() {
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscViewerDestroy(&_viewer);
}
void PetscHDF5::PushGroup(const std::string& group_name) {
    PetscViewerHDF5PushGroup(_viewer, group_name.c_str());
}
void PetscHDF5::PopGroup() {
    PetscViewerHDF5PopGroup(_viewer);
}
bool PetscHDF5::HasGroup(const std::string& group_name) const {
    PetscBool has;
    PetscErrorCode ierr;
    ierr = PetscViewerHDF5HasGroup(_viewer, group_name.c_str(), &has);
    return has;
}

void PetscHDF5::ReadAttribute(const std::string& attr_name, complex* value) {
    PetscErrorCode ierr;
    ierr = PetscViewerHDF5ReadAttribute(_viewer, "", attr_name.c_str(), PETSC_COMPLEX, value, value);
}
void PetscHDF5::ReadAttribute(const Vector object, const std::string& attr_name, complex* value) {
    PetscErrorCode ierr;
    auto petscVec = std::dynamic_pointer_cast<PetscVector>(object);
    ierr = PetscViewerHDF5ReadObjectAttribute(_viewer, (PetscObject)petscVec->_petsc_vec, attr_name.c_str(), PETSC_COMPLEX, value, value);
}

void PetscHDF5::ReadAttribute(const std::string& attr_name, int* value) {
    PetscErrorCode ierr;
    ierr = PetscViewerHDF5ReadAttribute(_viewer, "", attr_name.c_str(), PETSC_INT, value, value);
}
void PetscHDF5::ReadAttribute(const Vector object, const std::string& attr_name, int* value) {
    PetscErrorCode ierr;
    auto petscVec = std::dynamic_pointer_cast<PetscVector>(object);
    ierr = PetscViewerHDF5ReadObjectAttribute(_viewer, (PetscObject)petscVec->_petsc_vec, attr_name.c_str(), PETSC_INT, value, value);
}

void PetscHDF5::ReadAttribute(const std::string& attr_name, double* value) {
    PetscErrorCode ierr;
    ierr = PetscViewerHDF5ReadAttribute(_viewer, "", attr_name.c_str(), PETSC_DOUBLE, value, value);
}
void PetscHDF5::ReadAttribute(const Vector object, const std::string& attr_name, double* value) {
    PetscErrorCode ierr;
    auto petscVec = std::dynamic_pointer_cast<PetscVector>(object);
    ierr = PetscViewerHDF5ReadObjectAttribute(_viewer, (PetscObject)petscVec->_petsc_vec, attr_name.c_str(), PETSC_DOUBLE, value, value);
}



void PetscHDF5::WriteAttribute(const Vector object, const std::string& attr_name, const complex value) {
    auto petscVec = std::dynamic_pointer_cast<PetscVector>(object);
    PetscViewerHDF5WriteObjectAttribute(_viewer, (PetscObject)petscVec->_petsc_vec, attr_name.c_str(), PETSC_COMPLEX, &value);
}
void PetscHDF5::WriteAttribute(const std::string& attr_name, const complex value) {
    PetscViewerHDF5WriteAttribute(_viewer, "", attr_name.c_str(), PETSC_COMPLEX, &value);
}

void PetscHDF5::WriteAttribute(const Vector object, const std::string& attr_name, const double value) {
    auto petscVec = std::dynamic_pointer_cast<PetscVector>(object);
    PetscViewerHDF5WriteObjectAttribute(_viewer, (PetscObject)petscVec->_petsc_vec, attr_name.c_str(), PETSC_DOUBLE, &value);
}
void PetscHDF5::WriteAttribute(const std::string& attr_name, const double value) {
    PetscViewerHDF5WriteAttribute(_viewer, "", attr_name.c_str(), PETSC_DOUBLE, &value);
}
void PetscHDF5::WriteAttribute(const Vector object, const std::string& attr_name, const int value) {
    auto petscVec = std::dynamic_pointer_cast<PetscVector>(object);
    PetscViewerHDF5WriteObjectAttribute(_viewer, (PetscObject)petscVec->_petsc_vec, attr_name.c_str(), PETSC_INT, &value);
}
void PetscHDF5::WriteAttribute(const std::string& attr_name, const int value) {
    PetscViewerHDF5WriteAttribute(_viewer, "", attr_name.c_str(), PETSC_INT, &value);
}
bool PetscHDF5::HasAttribute(const std::string& attr_name) const {
    PetscBool has;
    PetscErrorCode ierr;
    ierr = PetscViewerHDF5HasAttribute(_viewer, NULL, attr_name.c_str(), &has);
    return has;
}
void PetscHDF5::WriteVector(const std::string& obj_name, const Vector value) {
    auto petscVec = std::dynamic_pointer_cast<PetscVector>(value);
    
    PetscObjectSetName((PetscObject)petscVec->_petsc_vec, obj_name.c_str());
    VecView(petscVec->_petsc_vec, _viewer);
}
void PetscHDF5::ReadVector(const std::string& obj_name, Vector value) {
    auto petscVec = std::dynamic_pointer_cast<PetscVector>(value);
    
    PetscObjectSetName((PetscObject)petscVec->_petsc_vec, obj_name.c_str());
    VecLoad(petscVec->_petsc_vec, _viewer);
}
bool PetscHDF5::HasVector(const std::string& obj_name) const {
    PetscBool has;
    PetscErrorCode ierr;
    ierr = PetscViewerHDF5HasDataset(_viewer, obj_name.c_str(), &has);

    // Vec b;
    // ierr = VecCreate(PETSC_COMM_WORLD,&b);PETSCASSERT(ierr);
    // ierr = PetscViewerHDF5HasObject(_viewer, (PetscObject)b, &has);
    // VecDestroy(&b);
    return has;
}
