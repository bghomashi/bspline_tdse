#pragma once

#include "matrix.h"
#include <vector>
#include <iostream>

class DenseMatrix : public Matrix {
    std::vector<complex> _mat;
    int _rows, _cols;
public:
    DenseMatrix() : _rows(0), _cols(0) {}
    DenseMatrix(int rows, int cols) : _rows(rows), _cols(cols), _mat(rows*cols) {}

    complex Get(int row, int col) const{
        return _mat[col + row*_cols];
    }
    void Set(int row, int col, complex val) {
        _mat[col + row*_cols] = val;
    }
    void Add(int row, int col, complex val) {
        _mat[col + row*_cols] += val;
    }

    void Scale(complex factor) {
        for (auto& e : _mat)
            e *= factor;
    }
    void Add(const Matrix& x, complex factor = 1.) {
        const DenseMatrix& dx = (const DenseMatrix&)x;
        for (int r : range(0,_rows)) {
            for (int c : range(0,_cols)) {
                _mat[c + r*_cols] += factor*dx._mat[c + r*_cols];
            }
        }
    }
    void Mult(const Vector& in, Vector& out) {

    }
    void CopyFrom(const Matrix& x) {
        const DenseMatrix& dx = (const DenseMatrix&)x;
        for (int r : range(0,_rows)) {
            for (int c : range(0,_cols)) {
                _mat[c + r*_cols] = dx._mat[c + r*_cols];
            }
        }
    }

    void print() const {
        for (int r : range(0,_rows)) {
            for (int c : range(0,_cols)) {
                std::cout << _mat[c + r*_cols] << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }
};