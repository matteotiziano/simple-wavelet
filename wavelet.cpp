#include "vector.h"
#include "wavelet.h"

using namespace Vector;

//=============================================================================
// Decomposition1D
//=============================================================================

template<typename T>
Decomposition1D<T>::Decomposition1D(size_t num_levels) {
  // Internal representation is structured as
  // {d1, d2, ... dn, an}
  this->resize(num_levels + 1);
}

template<typename T>
Decomposition1D<T>::~Decomposition1D() {}

template<typename T>
size_t Decomposition1D<T>::NumLevels() const {
  return this->size() - 1;
}

template<typename T>
const std::vector<T>& Decomposition1D<T>::GetAppcoef() const {
  return this->back();
}

template<typename T>
const std::vector<T>& Decomposition1D<T>::GetDetcoef(size_t level) const {
  if(level >= NumLevels()) {
    throw std::invalid_argument("Decomposition1D::GetDetcoef");
  }
  return this->at(level);
}

template<typename T>
void Decomposition1D<T>::SetDetcoef(const std::vector<T>& d, size_t level) {
  if(level >= NumLevels()) {
    throw std::invalid_argument("Decomposition1D::SetDetcoef");
  }
  this->operator[](level) = d;
}

template<typename T>
void Decomposition1D<T>::SetAppcoef(const std::vector<T>& a) {
  this->back() = a;
}

//=============================================================================
// Decomposition2D
//=============================================================================

template<typename T>
Decomposition2D<T>::Decomposition2D(size_t num_levels) {
  // Internal representation is structured as
  // {d1, v1, h1, d2, v2, h2, ... dn, vn, hn, an}
  this->resize(3 * num_levels + 1);
}

template<typename T>
Decomposition2D<T>::~Decomposition2D() {}

template<typename T>
size_t Decomposition2D<T>::NumLevels() const {
  return (this->size() - 1) / 3;
}

template<typename T>
void Decomposition2D<T>::SetAppcoef(const Matrix<T>& a) {
  this->back() = a;
}

template<typename T>
const Matrix<T>& Decomposition2D<T>::GetAppcoef() const {
  return this->back();
}

template<typename T>
const Matrix<T>& Decomposition2D<T>::GetDetcoef(WaveledSubbdand subband, size_t n) const {
  if (subband != kHorizontalSubband && subband != kVerticalSubband && subband != kDiagonalSubband) {
    throw std::invalid_argument("Decomposition2D::GetDetcoef");
  }
  if (n >= NumLevels()) {
    throw std::invalid_argument("Decomposition2D::GetDetcoef ");
  }
  // Convert subband to int for indexing
  const size_t idx = static_cast<size_t>(subband);
  return this->at(3 * n + idx);
}

template<typename T>
void Decomposition2D<T>::SetDetcoef(const Matrix<T>& d, WaveledSubbdand subband, size_t n) {
  if (subband != kHorizontalSubband && subband != kVerticalSubband && subband != kDiagonalSubband) {
    throw std::invalid_argument("Decomposition2D::SetDetcoef");
  }
  if (n >= NumLevels()) {
    throw std::invalid_argument("Decomposition2D::SetDetcoef ");
  }
  const size_t idx = static_cast<size_t>(subband);
  this->at(3 * n + idx) = d;
}

//=============================================================================
// Wavelet c-tor and d-tor
//=============================================================================

template<typename T>
Wavelet<T>::Wavelet(const std::vector<T>& Lo_D, const std::vector<T>& Hi_D,
                    const std::vector<T>& Lo_R, const std::vector<T>& Hi_R) {
  Lo_D_ = Lo_D;
  Hi_D_ = Hi_D;
  Lo_R_ = Lo_R;
  Hi_R_ = Hi_R;
  half_length_filter_ = Lo_D_.size() / 2;
}

template<typename T>
Wavelet<T>::~Wavelet()
{ }

//=============================================================================
// 1D Discrete Wavelet Transform (DWT)
//=============================================================================

template<typename T>
Decomposition1D<T> Wavelet<T>::Wavedec(const std::vector<T>& x, size_t num_levels) const {
  if (x.size() < (1 << num_levels)) {
    throw std::invalid_argument("Wavelet::Wavedec");
  }
  Decomposition1D<T> decomposition(num_levels);
  std::vector<T> a = x;
  for (size_t k = 0; k < num_levels; k++) {
    std::vector<T> ak, dk;
    Dwt(a, ak, dk);
    // Save approximation coefficients for next iteration
    a = ak;
    // Store detail coefficients of current level k to Decomposition object
    decomposition.SetDetcoef(dk, k);
  }
  // Store final approximation coefficients
  decomposition.SetAppcoef(a);
  return decomposition;
}

template<typename T>
std::vector<T> Wavelet<T>::Waverec(const Decomposition1D<T>& decomposition, size_t length_rec) const {
  return Appcoef(decomposition, 0, length_rec);
}

template<typename T>
void Wavelet<T>::Dwt(const std::vector<T>& x, std::vector<T>& a, std::vector<T>& d) const {
  a = Convdown(x, Lo_D_);
  d = Convdown(x, Hi_D_);
}

template<typename T>
std::vector<T> Wavelet<T>::Idwt(const std::vector<T>& a, const std::vector<T>& d, size_t l) const {
  return Upconv(a, Lo_R_, l) + Upconv(d, Hi_R_, l);
}

//=============================================================================
// 2D Discrete Wavelet Transform (DWT)
//=============================================================================

template<typename T>
Decomposition2D<T> Wavelet<T>::Wavedec(const Matrix<T>& x, size_t num_levels) const {
  if (x.GetNumRows() < (1 << num_levels) || x.GetNumColumns() < (1 << num_levels)) {
    throw std::invalid_argument("Wavelet::Wavedec");
  }
  Decomposition2D<T> decomposition(num_levels);
  Matrix<T> a = x;
  for (size_t k = 0; k < num_levels; k++) {
    Matrix<T> ak, hk, vk, dk;
    Dwt(a, ak, hk, vk, dk);
    // Save approximation coefficients for next iteration
    a = ak;
    // Store subbands obtained at current level k to Decomposition object
    decomposition.SetDetcoef(dk, kDiagonalSubband, k);
    decomposition.SetDetcoef(vk, kVerticalSubband, k);
    decomposition.SetDetcoef(hk, kHorizontalSubband, k);
  }
  // Store final approximation
  decomposition.SetAppcoef(a);
  return decomposition;
}

template<typename T>
Matrix<T> Wavelet<T>::Waverec(const Decomposition2D<T>& decomposition, const MatrixSize& size_rec) const {
  // Call Appcoef with maximum level of reconstruction (0)
  return Appcoef(decomposition, 0, size_rec);
}

template<typename T>
void Wavelet<T>::Dwt(const Matrix<T>& x, Matrix<T>& a, Matrix<T>& h, Matrix<T>& v, Matrix<T>& d) const {
  Convdown(x, Lo_D_, a, h);
  Convdown(x, Hi_D_, v, d);
}

template<typename T>
Matrix<T> Wavelet<T>::Idwt(const Matrix<T>& a, const Matrix<T>& h, const Matrix<T>& v,
                           const Matrix<T>& d, const MatrixSize& s) const {
  // Return sum of upsampled subbands
  return Upconv(a, Lo_R_, Lo_R_, s) + Upconv(h, Hi_R_, Lo_R_, s) +
      Upconv(v, Lo_R_, Hi_R_, s) + Upconv(d, Hi_R_, Hi_R_, s);
}

//=============================================================================
// Utils
//=============================================================================

template<typename T>
std::vector<T> Wavelet<T>::Linrec(const std::vector<T>& a, size_t num_levels) const {
  Decomposition1D<T> c(num_levels);
  const size_t length = a.size();
  for (size_t k = 0; k < num_levels; k++) {
    // Set k-th level detail coefficients to zero
    std::vector<T> dk(length * (1 << (num_levels - k - 1)), static_cast<T>(0));
    c.SetDetcoef(dk, k);
  }
  // Set approx coefficients equal to the given input
  c.SetAppcoef(a);
  // Call Appcoef with length of reconstruction equal to N*2^n
  return Appcoef(c, 0, length * (1 << num_levels));
}

template<typename T>
Matrix<T> Wavelet<T>::Linrec(const Matrix<T>& a, size_t num_levels) const {
  Decomposition2D<T> c(num_levels);
  const size_t num_rows = a.GetNumRows();
  const size_t num_cols = a.GetNumColumns();
  for (size_t k = 0; k < num_levels; k++) {
    // Set all subbands at k-th level to zero
    Matrix<T> tmp(num_rows * (1 << (num_levels - k - 1)), num_cols * (1 << (num_levels - k - 1)), static_cast<T>(0));
    c.SetDetcoef(tmp, kDiagonalSubband, k);
    c.SetDetcoef(tmp, kVerticalSubband, k);
    c.SetDetcoef(tmp, kHorizontalSubband, k);
  }
  // Set approx coefficients equal to given input
  c.SetAppcoef(a);
  // Call Appcoef with size_rec equal to to M*2^n-by-N*2^n
  const MatrixSize size_rec(num_rows * (1 << num_levels), num_cols * (1 << num_levels));
  return Appcoef(c, 0, size_rec);
}

//===========================================================================
// Utils for 1D DWT
//===========================================================================

template<typename T>
std::vector<T> Wavelet<T>::Appcoef(const Decomposition1D<T>& decomposition, size_t level) const {
  if (level >= decomposition.NumLevels()) {
    throw std::invalid_argument("Wavelet::Appcoef");
  }
  const size_t length_rec = (level == 0) ? 2 * decomposition.GetDetcoef(0).size() : decomposition.GetDetcoef(level - 1).size();
  return Appcoef(decomposition, level, length_rec);
}

template<typename T>
std::vector<T> Wavelet<T>::Appcoef(const Decomposition1D<T>& decomposition, size_t level, size_t length_rec) const {
  std::vector<T> a = decomposition.GetAppcoef();
  for (size_t k = decomposition.NumLevels(); k != level; k--) {
    // Get size from next upper level, if current level is the final then
    // resort to the final lenght argument sz
    const size_t length_rec_k = k != 1 ? decomposition.GetDetcoef(k - 2).size() : length_rec;
    // Reconstruct signal of size s from kth level approx and detail
    a = Idwt(a, decomposition.GetDetcoef(k - 1), length_rec_k);
  }
  return a;
}

template<typename T>
std::vector<T> Wavelet<T>::Convdown(const std::vector<T>& x, const std::vector<T>& f) const {
  const size_t last = static_cast<size_t>(std::ceil(static_cast<double>(x.size()) / 2.0) * 2.0);
  std::vector<T> z = Wextend(x, half_length_filter_);
  z = Conv(z, f, kValidConvolution);
  std::vector<T> y(last / 2);
  for (size_t k = 1; k <= last; k += 2) {
    // Since k is always an odd integer, k/2 === floor(k/2)
    y[k / 2] = z[k];
  }
  return y;
}

template<typename T>
std::vector<T> Wavelet<T>::Upconv(const std::vector<T>& x, const std::vector<T>& f, size_t s) const {
  const size_t lf = f.size();
  std::vector<T> y = Dyadup(x);
  y = Wextend(y, lf / 2);
  y = Conv(y, f, kFullConvolution);
  y.erase(y.begin(), y.begin() + lf - 1);
  y.erase(y.begin() + s, y.end());
  return y;
}

//===========================================================================
// Utils for 2D DWT
//===========================================================================

template<typename T>
Matrix<T> Wavelet<T>::Appcoef(const Decomposition2D<T>& decomposition, size_t level) const {
  if (level >= decomposition.NumLevels()) {
    throw std::invalid_argument("Wavelet::Appcoef2");
  }
  MatrixSize size_rec;
  if(level == 0) {
    // Level zero contains the biggest size of the decomposition, so we set the
    // reconstruction size equal to twice the size at level zero.
    const MatrixSize size_level_0 = decomposition.GetDetcoef(kHorizontalSubband, 0).Size();
    size_rec = MatrixSize(size_level_0.GetNumRows() * 2, size_level_0.GetNumColums() * 2);
  }
  else {
    // If level is not zero, get size of (any) subband at the next highest
    // level. We take here the horizontal subband.
    size_rec = decomposition.GetDetcoef(kHorizontalSubband, level - 1).Size();
  }
  return Appcoef(decomposition, level, size_rec);
}

template<typename T>
Matrix<T> Wavelet<T>::Appcoef(const Decomposition2D<T>& decomposition, size_t level,
                              const MatrixSize& size_rec) const {
  Matrix<T> a = decomposition.GetAppcoef();
  for (size_t k = decomposition.NumLevels(); k != level; k--) {
    // Select size of reconstruction from any subband at the next highest level
    // We arbitrarily use the horizontal subband. If there is no next highest
    // level, then we use the size_rec argument.
    const MatrixSize s = (k != 1) ? decomposition.GetDetcoef(kHorizontalSubband, k - 2).Size() : size_rec;
    const Matrix<T>& hk = decomposition.GetDetcoef(kHorizontalSubband, k - 1);
    const Matrix<T>& vk = decomposition.GetDetcoef(kVerticalSubband, k - 1);
    const Matrix<T>& dk = decomposition.GetDetcoef(kDiagonalSubband, k - 1);
    a = Idwt(a, hk, vk, dk, s);
  }
  return a;
}

template<typename T>
void Wavelet<T>::Convdown(const Matrix<T>& x, const std::vector<T>& f, Matrix<T>& xL, Matrix<T>& xH) const {
  // Compute size of subbands
  const size_t num_rows_subband = static_cast<size_t>(std::ceil(static_cast<double>(x.GetNumRows()) / 2.0));
  const size_t num_cols_subband = static_cast<size_t>(std::ceil(static_cast<double>(x.GetNumColumns()) / 2.0));
  // Convolve and decimate each row with the chosen filter
  Matrix<T> tmp(x.GetNumRows(), num_cols_subband);
  for (size_t i = 0; i < x.GetNumRows(); i++) {
    const std::vector<T> row = Convdown(x.GetRow(i), f);
    tmp.SetRow(row, i);
  }
  // Compute high-pass subband of pre-filtered tmp matrix
  xH = Matrix<T>(num_rows_subband, num_cols_subband);
  for (size_t j = 0; j < tmp.GetNumColumns(); j++) {
    const std::vector<T> column = Convdown(tmp.GetColumn(j), Hi_D_);
    xH.SetColumn(column, j);
  }
  // Compute low-pass subband of pre-filtered tmp matrix
  xL = Matrix<T>(num_rows_subband, num_cols_subband);
  for (size_t j = 0; j < tmp.GetNumColumns(); j++) {
    const std::vector<T> column = Convdown(tmp.GetColumn(j), Lo_D_);
    xL.SetColumn(column, j);
  }
}

template<typename T>
Matrix<T> Wavelet<T>::Upconv(const Matrix<T>& x, const std::vector<T>& f1, const std::vector<T>& f2, const MatrixSize& s) const {
  // Upsample each column
  Matrix<T> tmp(s[0], x.GetNumColumns());
  for (size_t j = 0; j < x.GetNumColumns(); j++) {
    tmp.SetColumn(Upconv(x.GetColumn(j), f1, s.GetNumRows()), j);
  }
  // Upsample each row of tmp Matrix
  Matrix<T> y(s.GetNumRows(), s.GetNumColums());
  for (size_t i = 0; i < y.GetNumRows(); i++) {
    y.SetRow(Upconv(tmp.GetRow(i), f2, s.GetNumColums()), i);
  }
  return y;
}

//=============================================================================
// Utils
//=============================================================================

template<typename T>
std::vector<T> Wavelet<T>::Dyadup(const std::vector<T>& x) const {
  // Return a vector having double length
  std::vector<T> y;
  y.reserve(2 * x.size());
  // The output vector interleaves the elements of x with zero
  for (size_t k = 0; k < x.size(); k++) {
    y.push_back(x[k]);
    y.push_back(0);
  }
  return y;
}

template<typename T>
std::vector<T> Wavelet<T>::Wextend(const std::vector<T>& x, size_t lenEXT) const {
  std::vector<T> temp = x;
  // Add extra sample if signal is odd
  if (x.size() % 2 == 1) {
    temp.push_back(temp.back());
  }
  // Handle cases when x is shorter than lenEXT
  const size_t rep = lenEXT / temp.size(); // (size_t)std::floor(lenEXT/temp.size());
  const size_t len = lenEXT % temp.size();
  std::vector<T> y;
  y.reserve(2 * (len + rep * temp.size()) + temp.size());
  // Copy last len elements at the beginning
  y.insert(y.begin(), temp.end() - len, temp.end());
  for (size_t k = 0; k < 2 * rep + 1; k++) {
    y.insert(y.end(), temp.begin(), temp.end());
  }
  // Copy first len elements at the end
  y.insert(y.end(), temp.begin(), temp.begin() + len);
  return y;
}

// Explicit template instantiation
template class Decomposition1D<double>;
template class Decomposition2D<double>;
template class Wavelet<double>;
