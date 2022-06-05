
#pragma once

#include <complex>
#include <functional>
#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

class IMatrix;
class IVector;
class IHDF5;
class IGMRESSolver;
class IASCII;

typedef std::complex<double> complex;
typedef std::function<complex(int, int)> FuncOfRowCol;
typedef std::shared_ptr<IMatrix> Matrix;
typedef std::shared_ptr<IVector> Vector;
typedef std::shared_ptr<IHDF5> HDF5;
typedef std::shared_ptr<IASCII> ASCII;
typedef std::shared_ptr<IGMRESSolver> GMRESSolver;

const double c = 137.036;
const double Pi = 3.14159265358979323846;
const double LnmToEnergy = 45.5633346744;


enum DimIndex {
    X = 0,
    Y,
    Z,

    NUM
};

template <typename T>
std::vector<T> range(T start, T end, T step = 1) {
    std::vector<T> r; r.reserve((end - start)/step);
    for (T i = start; i < end; i+=step)
        r.push_back(i);
    return r;
}

class IVector {
protected:
    int _len;
public:
    inline int Length() const {
        return _len;
    }

    virtual complex Get(int index) const = 0;
    virtual void Get(std::vector<complex>& out) = 0;
    virtual void Set(int index, complex value) = 0;
    virtual void Scale(complex a) = 0;
    virtual void Duplicate(const Vector& o) = 0;
    virtual void Copy(const Vector& o) = 0;
    virtual void Zero() = 0;
    virtual void Concatenate(const std::vector<Vector>& vecs) = 0;
    virtual void CopyTo(std::vector<complex>& values) = 0;
    virtual void Transform(Vector& out, std::function<std::vector<complex>(const std::vector<complex>&)> f) = 0;

    virtual Vector GetSubVector(int start, int end) = 0;
    virtual void RestoreSubVector(Vector sub) = 0; 
    virtual void AssembleBegin() {};
    virtual void AssembleEnd() {};
};

class IMatrix {
protected:
    int _rows, _cols;
public:
    int Rows() const {
        return _rows;
    }
    int Cols() const {
        return _cols;
    }

    virtual complex Get(int row, int col) const = 0;
    virtual void Set(int row, int col, complex val) = 0;
    virtual void Add(int row, int col, complex val) = 0;
    virtual void Mult(const Vector in, Vector out) = 0;
    virtual void Scale(complex factor) = 0;
    virtual void Duplicate(const Matrix o) = 0;               // allocate AND copy values
    virtual void Zero() = 0;
    virtual void Copy(const Matrix o) = 0;                    // just copy values

    virtual void AssembleBegin() {};
    virtual void AssembleEnd() {};
    
    virtual bool IsSymmetric(double tol = 1e-17) const = 0;
    virtual bool IsAntiSymmetric(double tol = 1e-17) const = 0;

    virtual void FillBandedBlock(int bandwidth, int blocksize, int blockRow, int blockCol, FuncOfRowCol element) = 0;
    virtual void FillBandedBlock(int bandwidth, int blocksize, FuncOfRowCol element) = 0;
    virtual void FillBandedBlock(int bandwidth, int blockSize, int blockColOffset, FuncOfRowCol element) = 0;
};

class IASCII {
protected:
    std::fstream _file;
public:
    IASCII();
    bool Open(const std::string& filename, char mode);
    void Close();

    virtual void Flush();
    virtual void Write(const std::string& text);
    virtual std::string ReadLine();
};

class IHDF5 {
public:
    virtual void PushGroup(const std::string& group_name) = 0;
    virtual void PopGroup() = 0;
    virtual bool HasGroup(const std::string& group_name) const = 0;

    virtual void ReadAttribute(const std::string& attr_name, complex* value) = 0;
    virtual void ReadAttribute(const std::string& attr_name, double* value) = 0;
    virtual void ReadAttribute(const std::string& attr_name, int* value) = 0;
    virtual void ReadAttribute(const Vector object, const std::string& attr_name, complex* value) = 0;
    virtual void ReadAttribute(const Vector object, const std::string& attr_name, double* value) = 0;
    virtual void ReadAttribute(const Vector object, const std::string& attr_name, int* value) = 0;
    virtual void WriteAttribute(const Vector object, const std::string& attr_name, const complex value) = 0;
    virtual void WriteAttribute(const Vector object, const std::string& attr_name, const double value) = 0;
    virtual void WriteAttribute(const Vector object, const std::string& attr_name, const int value) = 0;
    virtual void WriteAttribute(const std::string& attr_name, const complex value) = 0;
    virtual void WriteAttribute(const std::string& attr_name, const double value) = 0;
    virtual void WriteAttribute(const std::string& attr_name, const int value) = 0;
    virtual bool HasAttribute(const std::string& attr_name) const = 0;
    virtual void WriteVector(const std::string& obj_name, const Vector value) = 0;
    virtual void ReadVector(const std::string& obj_name, Vector value) = 0;
    virtual bool HasVector(const std::string& obj_name) const = 0;
};

class IGMRESSolver {
public:
    virtual bool Solve(const Matrix A, const Vector b, Vector x) = 0;
    virtual void SetBlockedPC(int blocks) = 0;
};



class MathLib {
public:
    virtual bool Startup(int argc, char **args) = 0;
    virtual void Shutdown() = 0;

    virtual Vector CreateVector(int N) = 0;
    virtual void DestroyVector(Vector& m) = 0;

    virtual Matrix CreateMatrix(int rows, int cols, int numBands) = 0;
    virtual void DestroyMatrix(Matrix& m) = 0;

    virtual GMRESSolver CreateGMRESSolver(int restart_iter = 500, int max_iter = 10000) = 0;
    virtual void DestroyGMRESSolver(GMRESSolver& m) = 0;

    virtual void Mult(const Matrix M, const Vector in, Vector out) = 0;
    virtual void Dot(const Vector a, const Vector b, complex& value) = 0;
    virtual void AYPX(Matrix Y, complex a, const Matrix X) = 0;
    virtual void AXPY(Matrix Y, complex a, const Matrix X) = 0;
    virtual void AYPX(Vector Y, complex a, const Vector X) = 0;
    virtual void AXPY(Vector Y, complex a, const Vector X) = 0;

    virtual void Eigen( const Matrix A, const Matrix S, 
                int numVectors, double tol,
                std::vector<complex>& values, 
                std::vector<Vector>& vectors) = 0;

    virtual HDF5 OpenHDF5(const std::string& filename, char mode) = 0;
    virtual void CloseHDF5(HDF5& file) = 0;

    virtual ASCII OpenASCII(const std::string& filename, char mode) = 0;
    virtual void CloseASCII(ASCII& file) = 0;

    virtual void ParallelPrintf(const std::string& text) = 0;
};



