#ifndef MATRIX_H_
#define MATRIX_H_

#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>

#include "vector.h"

//=============================================================================
// MatrixSize
//=============================================================================

struct MatrixSize : public std::array<size_t, 2> {
  MatrixSize() {}
  MatrixSize(size_t num_rows, size_t num_cols) {
    this->operator[](0) = num_rows;
    this->operator[](1) = num_cols;
  }
  size_t GetNumRows() const { return this->at(0); }
  size_t GetNumColums() const { return this->at(1); }
};

//=============================================================================
// Matrix
// Simple data type that allows to represent a 2D matrix in column-major
// order. This class also offers several operators to manipulate the data.
//=============================================================================

template<typename T>
class Matrix {

public:

  //===========================================================================
  // Matrix c-tors and d-tors
  //===========================================================================

  Matrix(size_t num_rows, size_t num_cols, T value = static_cast<T>(0));
  Matrix(size_t num_rows, size_t num_cols, const std::vector<T>& data);
  Matrix();
  ~Matrix();

  //===========================================================================
  // Basic accessors
  //===========================================================================

  size_t GetNumRows() const;
  size_t GetNumColumns() const;
  MatrixSize Size() const;
  size_t Numel() const;

  //===========================================================================
  // Indexing
  //===========================================================================

  std::array<size_t, 2> Index2Subscripts(size_t idx) const;
  size_t Subscripts2Index(size_t i, size_t j) const;

  T& operator()(size_t i, size_t j);
  const T& operator()(size_t i, size_t j) const;
  T& operator()(size_t idx);
  const T& operator()(size_t idx) const;

  //===========================================================================
  // Accessors
  //===========================================================================

  void SetRow(const std::vector<T>& row, size_t i);
  void SetColumn(const std::vector<T>& column, size_t j);
  void SetData(const std::vector<T>& data);

  std::vector<T> GetRow(size_t i) const;
  std::vector<T> GetColumn(size_t j) const;
  std::vector<T>& GetData();
  const std::vector<T>& GetData() const;

  void DeleteRow(size_t i);
  void DeleteColumn(size_t j);

  //===========================================================================
  // Utils
  //===========================================================================

  std::vector<T> Diag() const;

  void Fill(T value);
  void Rand();
  void Identity();

  void Transpose();
  void Reshape(size_t num_rows, size_t num_cols);

  void Round();
  void Floor();
  void Ceil();
  void Clamp(T min_value, T max_value);

  T Max() const;
  T Min() const;

  //===========================================================================
  // Binary Matrix-Matrix operations
  //===========================================================================

  Matrix<T> operator*(const Matrix<T>& rhs) const;
  Matrix<T> operator+(Matrix<T> rhs) const;
  Matrix<T> operator-(Matrix<T> rhs) const;

  Matrix<T> ElementProduct(Matrix<T> rhs) const;
  Matrix<T> ElementQuotient(Matrix<T> rhs) const;

  //===========================================================================
  // Binary Scalar-Matrix operations
  //===========================================================================

  Matrix<T> operator+(T rhs) const;
  Matrix<T> operator-(T rhs) const;
  Matrix<T> operator*(T rhs) const;
  Matrix<T> operator/(T rhs) const;
  Matrix<T> operator^(T rhs) const;

  //===========================================================================
  // Unary operations
  //===========================================================================

  Matrix<T> operator-() const;

private:

  //===========================================================================
  // Member variables
  //===========================================================================

  size_t num_rows_;
  size_t num_cols_;
  std::vector<T> data_;

}; // class Matrix<T>

//=============================================================================
// Overloading operator<< for ostream operations
//=============================================================================

// Overloaded operator<< for Matrix<T>
template<typename T>
inline std::ostream& operator<<(std::ostream& out, const Matrix<T>& matrix) {
  for (size_t i = 0, endI = matrix.GetNumRows(); i < endI; i++) {
    out << matrix.GetRow(i) << std::endl;
  }
  return out;
}

#endif // MATRIX_H_
