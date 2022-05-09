#pragma once

#include "maths/maths.h"
#include "utility/logger.h"
#include "utility/profiler.h"

#include <cassert>
#include <petsc.h>
#include <slepc.h>

#define PETSCASSERT(x) assert(x == 0)

class Petsc;
class PetscVector;
class PetscMatrix;
class EPSSolver;
class KSPSolver;

class PetscVector : public IVector {
    friend Petsc;
    friend EPSSolver;
    friend KSPSolver;
    friend void MatMult(const Matrix& M, const Vector& in, Vector& out);
    friend void MatAXPY(Matrix& X, complex a, const Matrix& Y);

    bool _dirty;
    std::vector<complex> _local_copy;
public:
    typedef std::shared_ptr<PetscVector> Ptr_t;

    VecScatter _petsc_ctx;
    Vec _petsc_sca_vec;
    Vec _petsc_vec;
    IS _petsc_is;

    PetscVector();
    PetscVector(int length);
    ~PetscVector();

    void AssembleBegin();
    void AssembleEnd();

    complex Get(int index) const;
    void Get(std::vector<complex>& out);
    void Set(int index, complex value);
    void Scale(complex a);
    void Duplicate(const Vector& o);
    void Copy(const Vector& o);
    void Zero();
    void Concatenate(const std::vector<Vector>& vecs);
    void CopyTo(std::vector<complex>& values); 
    void Transform(Vector& out, std::function<std::vector<complex>(const std::vector<complex>&)> f);
    Vector GetSubVector(int start, int end);
    void RestoreSubVector(Vector sub);

    void CreateScatter();
    void Scatter();
    void ScatterGetArray(complex** ptr);
    void ScatterRestoreArray(complex** ptr);
};

class PetscMatrix : public IMatrix {
    friend Petsc;
    friend EPSSolver;
    friend KSPSolver;
    friend void MatMult(const Matrix& M, const Vector& in, Vector& out);
    friend void MatAXPY(Matrix& X, complex a, const Matrix& Y);

public:
    typedef std::shared_ptr<PetscMatrix> Ptr_t;

    Mat _petsc_mat;
    int _row_start, _row_end;

    PetscMatrix();
    PetscMatrix(int rows, int cols, int numbands);
    PetscMatrix(const PetscMatrix& o);
    ~PetscMatrix();

    void AssembleBegin();
    void AssembleEnd();


    complex Get(int row, int col) const;
    void Set(int row, int col, complex val);
    void Add(int row, int col, complex val);
    void Mult(const Vector in, Vector out);
    void Scale(complex factor);
    void Duplicate(const Matrix o);               // allocate AND copy values
    void Zero();
    void Copy(const Matrix o);      
    

    void FillBandedBlock(int bandwidth, int blocksize, int blockRow, int blockCol, FuncOfRowCol element);
    void FillBandedBlock(int bandwidth, int blocksize, FuncOfRowCol element);
    void FillBandedBlock(int bandwidth, int blockSize, int blockColOffset, FuncOfRowCol element);


    void Set(const std::vector<int>& rows, const std::vector<int>& cols, const complex* value);
    
    bool IsSymmetric(double tol) const;
    bool IsAntiSymmetric(double tol) const;

    void HermitianTranspose(Matrix& transpose);
    void Transpose(Matrix& transpose);

    // ------- petsc specific --------
    int RowStart() const;
    int RowEnd() const;
        
};

class PetscASCII : public IASCII {
    
public:
    PetscASCII();
    PetscASCII(const std::string& filename, char mode);
    ~PetscASCII();
 
    void Write(const std::string& text);
    std::string ReadLine();
    void Flush();
};

class PetscHDF5 : public IHDF5 {
    PetscViewer _viewer;
public:
    PetscHDF5();
    PetscHDF5(const std::string& filename, char mode);
    ~PetscHDF5();

    void PushGroup(const std::string& group_name);
    void PopGroup();

    void WriteAttribute(const Vector object, const std::string& attr_name, const complex value);
    void WriteAttribute(const Vector object, const std::string& attr_name, const double value);
    void WriteAttribute(const Vector object, const std::string& attr_name, const int value);
    void WriteAttribute(const std::string& attr_name, const complex value);
    void WriteAttribute(const std::string& attr_name, const double value);
    void WriteAttribute(const std::string& attr_name, const int value);
    void WriteVector(const std::string& obj_name, const Vector value);
    void ReadVector(const std::string& obj_name, Vector value);
};

class PetscSolver : public IGMRESSolver {
    KSPConvergedReason _reason;
    const char *_strreason;
public:
    KSP _petsc_ksp;          /* linear solver context */
    PC _petsc_pc;            /* preconditioner context */

    PetscSolver(int restart_iter = 500, int max_iter = 10000);
    ~PetscSolver();

    void SetBlockedPC(int blocks);
    bool Solve(const Matrix A, const Vector b, Vector x);
};

class PetscLogger : public Logger {
public:
    void info(const std::string& text);
    void warn(const std::string& text);
    void critical(const std::string& text);
    void debug(const std::string& text);
    void set_logger_file(const std::string& log_file);
};

class PetscProfiler : public Profiler {
public:
    void Push(const std::string& name);
    void Pop(const std::string& name);
    void Print();
    bool PrintTo(const std::string& log_file);
};

class Petsc : public MathLib {
    // ksp
    PetscMPIInt _size, _rank;
public:
    bool Startup(int argc, char **args);
    void Shutdown();
    
    Vector CreateVector(int N);
    void DestroyVector(Vector& m);

    Matrix CreateMatrix(int rows, int cols, int numBands);
    void DestroyMatrix(Matrix& m);

    GMRESSolver CreateGMRESSolver(int restart_iter = 500, int max_iter = 10000);
    void DestroyGMRESSolver(GMRESSolver& m);

    HDF5 OpenHDF5(const std::string& filename, char mode);
    void CloseHDF5(HDF5& file);
    
    ASCII OpenASCII(const std::string& filename, char mode);
    void CloseASCII(ASCII& file);

    void Mult(const Matrix M, const Vector in, Vector out);
    void Dot(const Vector a, const Vector b, complex& value);
    void AYPX(Matrix Y, complex a, const Matrix X);
    void AXPY(Matrix Y, complex a, const Matrix X);
    void AYPX(Vector Y, complex a, const Vector X);
    void AXPY(Vector Y, complex a, const Vector X);

    void Eigen( const Matrix A, const Matrix S, 
                int numVectors, double tol,
                std::vector<complex>& values, 
                std::vector<Vector>& vectors);


    void ParallelPrintf(const std::string& text);


    // singleton - only one petsc
    static Petsc& get() {
        static Petsc sInstance;
        return sInstance;
    }
};
