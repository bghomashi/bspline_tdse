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
    PetscViewerDestroy(&_viewer);
    MPI_Barrier(PETSC_COMM_WORLD);
}
void PetscHDF5::PushGroup(const std::string& group_name) {
    PetscViewerHDF5PushGroup(_viewer, group_name.c_str());
}
void PetscHDF5::PopGroup() {
    PetscViewerHDF5PopGroup(_viewer);
}
// doesnt work??
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
// bool PetscVector::Load(const std::string& filename, const std::string& object_name, const std::string& group_name) {
//     PetscErrorCode ierr;
//     PetscViewer viewer;
    
//     ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"hdf5output.h5",FILE_MODE_READ,&viewer); CHKERRQ(ierr);
//     PetscViewerSetFromOptions(viewer);

//     PetscViewerHDF5SetBaseDimension2(viewer, PETSC_FALSE);
//     PetscObjectSetName((PetscObject) _petsc_vec, object_name.c_str());
//     VecLoad(_petsc_vec,viewer);
    
//     PetscViewerDestroy(&viewer);
//     MPI_Barrier(PETSC_COMM_WORLD);

//     return true;
// }

