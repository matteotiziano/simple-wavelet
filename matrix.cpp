#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <random>

#include "matrix.h"

using namespace Vector;

//=============================================================================
// Matrix c-tors and d-tor
//=============================================================================

template<typename T>
Matrix<T>::Matrix(size_t num_rows, size_t num_cols, T value) {
  num_rows_ = num_rows;
  num_cols_ = num_cols;
  data_.resize(num_rows_ * num_cols_, value);
}

template<typename T>
Matrix<T>::Matrix(size_t num_rows, size_t num_cols, const std::vector<T>& data) {
  if (num_rows * num_cols != data.size()) {
    throw std::invalid_argument("Matrix::Matrix");
  }
  num_rows_ = num_rows;
  num_cols_ = num_cols;
  data_ = data;
}

template<typename T>
Matrix<T>::Matrix() { }

template<typename T>
Matrix<T>::~Matrix() { }

//=============================================================================
// Basic accessors
//=============================================================================

template<typename T>
size_t Matrix<T>::GetNumRows() const {
  return num_rows_;
}

template<typename T>
size_t Matrix<T>::GetNumColumns() const {
  return num_cols_;
}

template<typename T>
MatrixSize Matrix<T>::Size() const {
  if (num_rows_ * num_cols_ != data_.size()) {
    throw std::logic_error("Matrix::Size");
  }
  // clang++ suggests braces around initialization of subobject
  return MatrixSize(num_rows_, num_cols_);
}

template<typename T>
size_t Matrix<T>::Numel() const {
  if (num_rows_ * num_cols_ != data_.size()) {
    throw std::logic_error("Matrix::Numel");
  }
  return data_.size();
}

//=============================================================================
// Indexing
//=============================================================================

template<typename T>
std::array<size_t, 2> Matrix<T>::Index2Subscripts(size_t idx) const {
  if (idx >= data_.size()) {
    throw std::out_of_range("Matrix::Index2Subscripts");
  }
  return {{idx % num_rows_, idx / num_rows_}};
}

template<typename T>
size_t Matrix<T>::Subscripts2Index(size_t i, size_t j) const {
  if (i >= num_rows_ || j >= num_cols_) {
    throw std::out_of_range("Matrix::Subscripts2Index");
  }
  return j * num_rows_ + i;
}

template<typename T>
T& Matrix<T>::operator()(size_t i, size_t j) {
  const size_t idx = Subscripts2Index(i, j);
  return operator()(idx);
}

template<typename T>
const T& Matrix<T>::operator()(size_t i, size_t j) const {
  const size_t idx = Subscripts2Index(i, j);
  return operator()(idx);
}

template<typename T>
T& Matrix<T>::operator()(size_t idx) {
  if (idx >= Numel()) {
    throw std::out_of_range("Matrix::operator()");
  }
  return data_[idx];
}

template<typename T>
const T& Matrix<T>::operator()(size_t idx) const {
  if (idx >= Numel()) {
    throw std::out_of_range("Matrix::operator()");
  }
  return data_.at(idx);
}

//=============================================================================
// Accessors
//=============================================================================

template<typename T>
void Matrix<T>::SetRow(const std::vector<T>& row, size_t i) {
  if (i >= num_rows_ || num_cols_ != row.size()) {
    throw std::invalid_argument("Matrix::SetRow");
  }
  for (size_t j = 0; j < num_cols_; j++) {
    operator()(i, j) = row.at(j);
  }
}

template<typename T>
void Matrix<T>::SetColumn(const std::vector<T>& column, size_t j) {
  if (j >= num_cols_ || num_rows_ != column.size()) {
    throw std::invalid_argument("Matrix::SetColumn");
  }
  const ptrdiff_t idx = static_cast<ptrdiff_t>(Subscripts2Index(0, j));
  std::copy(column.begin(), column.end(), data_.begin() + idx);
}

template<typename T>
void Matrix<T>::SetData(const std::vector<T>& data) {
  if (Numel() != data.size()) {
    throw std::invalid_argument("Matrix::SetData");
  }
  data_ = data;
}

template<typename T>
std::vector<T> Matrix<T>::GetRow(size_t i) const {
  std::vector<T> row;
  row.reserve(num_cols_);
  for (size_t j = 0; j < num_cols_; j++) {
    const T& value = operator()(i, j);
    row.push_back(value);
  }
  return row;
}

template<typename T>
std::vector<T> Matrix<T>::GetColumn(size_t j) const {
  const ptrdiff_t idx_start = static_cast<ptrdiff_t>(Subscripts2Index(0, j));
  const ptrdiff_t idx_end = idx_start + static_cast<ptrdiff_t>(num_rows_);
  return std::vector<T>(data_.begin() + idx_start, data_.begin() + idx_end);
}

template<typename T>
std::vector<T>& Matrix<T>::GetData() {
  return data_;
}

template<typename T>
const std::vector<T>& Matrix<T>::GetData() const {
  return data_;
}

template<typename T>
void Matrix<T>::DeleteRow(size_t i) {
  for (size_t j = 0; j < num_cols_; j++) {
    const ptrdiff_t idx = static_cast<ptrdiff_t>(Subscripts2Index(i, num_cols_ - j - 1));
    data_.erase(data_.begin() + idx);
  }
  num_rows_--;
  assert(Numel() == data_.size());
}

template<typename T>
void Matrix<T>::DeleteColumn(size_t j) {
  const ptrdiff_t idx_start = static_cast<ptrdiff_t>(Subscripts2Index(0, j));
  const ptrdiff_t idx_end = static_cast<ptrdiff_t>(Subscripts2Index(0, j + 1));
  data_.erase(data_.begin() + idx_start, data_.begin() + idx_end);
  num_cols_--;
  assert(Numel() == data_.size());
}

//=============================================================================
// Utils
//=============================================================================

template<typename T>
std::vector<T> Matrix<T>::Diag() const {
  std::vector<T> diag;
  const size_t diag_length = std::min(num_rows_, num_cols_);
  diag.reserve(diag_length);
  for (size_t k = 0; k < diag_length; k++) {
    const T& value = operator()(k, k);
    diag.push_back(value);
  }
  return diag;
}

template<typename T>
void Matrix<T>::Fill(T value) {
  std::fill(data_.begin(), data_.end(), value);
}

template<typename T>
void Matrix<T>::Rand(unsigned int seed) {
  std::mt19937 gen;
  gen.seed(seed);
  std::uniform_real_distribution<T> dis(0, 1);
  std::generate(data_.begin(), data_.end(), [&dis, &gen]() {
    return dis(gen);
  });
}

template<typename T>
void Matrix<T>::Identity() {
  if (num_rows_ != num_cols_) {
    throw std::logic_error("Matrix::Identity");
  }
  for (size_t j = 0; j < num_cols_; j++) {
    for (size_t i = 0; i < num_rows_; i++) {
      operator()(i, j) = (i == j)? static_cast<T>(1) : static_cast<T>(0);
    }
  }
}

template<typename T>
void Matrix<T>::Transpose() {
  // TODO: inplace tranposition
  std::vector<T> transposed_data;
  transposed_data.reserve(num_rows_ * num_cols_);
  for (size_t i = 0; i < num_rows_; i++) {
    for (size_t j = 0; j < num_cols_; j++) {
      const T& value = operator()(i, j);
      transposed_data.push_back(value);
    }
  }
  data_ = transposed_data;
  std::swap(num_rows_, num_cols_);
  assert(Numel() == data_.size());
}

template<typename T>
void Matrix<T>::Reshape(size_t num_rows, size_t num_cols) {
  if (num_rows * num_cols != Numel()) {
    throw std::invalid_argument("Matrix::Reshape");
  }
  num_rows_ = num_rows;
  num_cols_ = num_cols;
  assert(Numel() == data_.size());
}

template<typename T>
void Matrix<T>::Round() {
  Vector::Round(data_);
}

template<typename T>
void Matrix<T>::Floor() {
  Vector::Floor(data_);
}

template<typename T>
void Matrix<T>::Ceil() {
  Vector::Ceil(data_);
}

template<typename T>
void Matrix<T>::Clamp(T min_value, T max_value) {
  Vector::Clamp(data_, min_value, max_value);
}

template<typename T>
T Matrix<T>::Max() const {
  return Vector::Max(data_);
}

template<typename T>
T Matrix<T>::Min() const {
  return Vector::Min(data_);
}

//=============================================================================
// Binary Matrix-Matrix operations
//=============================================================================

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs) const {
  if (rhs.GetNumRows() != num_cols_) {
    throw std::invalid_argument("Matrix::operator*");
  }
  Matrix<T> res(num_rows_, rhs.GetNumColumns(), 0.0);
  for(size_t i = 0; i < num_rows_; i++) {
    for(size_t j = 0; j < rhs.GetNumColumns(); j++) {
      T& value = res(i, j);
      for(size_t k = 0; k < num_cols_; k++) {
        value += operator()(i, k) * rhs(k, j);
      }
    }
  }
  return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(Matrix<T> rhs) const {
  if (Size() != rhs.Size()) {
    throw std::invalid_argument("Matrix::operator+");
  }
  rhs.GetData() += data_;
  return rhs;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(Matrix<T> rhs) const {
  if (Size() != rhs.Size()) {
    throw std::invalid_argument("Matrix::operator-");
  }
  rhs.GetData() -= data_;
  return rhs;
}

template<typename T>
Matrix<T> Matrix<T>::ElementProduct(Matrix<T> rhs) const {
  if (Size() != rhs.Size()) {
    throw std::invalid_argument("Matrix::ElementProduct");
  }
  rhs.GetData() *= data_;
  return rhs;
}

template<typename T>
Matrix<T> Matrix<T>::ElementQuotient(Matrix<T> rhs) const {
  if (Size() != rhs.Size()) {
    throw std::invalid_argument("Matrix::ElementQuotient");
  }
  rhs.GetData() /= data_;
  return rhs;
}

//=============================================================================
// Binary Scalar-Matrix operations
//=============================================================================

template<typename T>
Matrix<T> Matrix<T>::operator+(T rhs) const {
  return Matrix<T>(num_rows_, num_cols_, data_ + rhs);
}

template<typename T>
Matrix<T> Matrix<T>::operator-(T rhs) const {
  return Matrix<T>(num_rows_, num_cols_, data_ - rhs);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(T rhs) const {
  return Matrix<T>(num_rows_, num_cols_, data_ * rhs);
}

template<typename T>
Matrix<T> Matrix<T>::operator/(T rhs) const {
  return Matrix<T>(num_rows_, num_cols_, data_ / rhs);
}

template<typename T>
Matrix<T> Matrix<T>::operator^(T rhs) const {
  return Matrix<T>(num_rows_, num_cols_, data_ ^ rhs);
}

//=============================================================================
// Unary operations
//=============================================================================

template<typename T>
Matrix<T> Matrix<T>::operator-() const {
  return Matrix<T>(num_rows_, num_cols_, -data_);
}

// Explicit template instantiation
template class Matrix<double>;
