#include <iostream>

#include "wavelet.h"

using namespace Vector;

int main(int argc, const char *argv[]) {

  // TODO: take parameters from argv

  // Initialise a random 2D matrix
  const size_t num_rows = 64;
  const size_t num_cols = 64;
  Matrix<double> matrix(num_rows, num_cols);
  matrix.Rand();

  // Definition of Haar filters
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  const std::vector<double> Lo_D = {inv_sqrt2, inv_sqrt2};
  const std::vector<double> Hi_D = {-inv_sqrt2, inv_sqrt2};
  const std::vector<double> Lo_R = {inv_sqrt2, inv_sqrt2};
  const std::vector<double> Hi_R = {inv_sqrt2, -inv_sqrt2};

  // Instantiate a Wavelet object with Haar filters
  const Wavelet<double> haar(Lo_D, Hi_D, Lo_R, Hi_R);

  // Number of decomposition levels
  const size_t num_levels = 4;

  // Example of 1D decomposition and reconstruction using the first column of
  // the above matrix as input
  std::vector<double> vector = matrix.GetColumn(0);
  const Decomposition1D<double> dec1D = haar.Wavedec(vector, num_levels);
  std::vector<double> vector_rec = haar.Waverec(dec1D, vector.size());
  std::cout << "Error of 1D reconstruction: " << Norm(vector - vector_rec) << std::endl;

  // Example of 2D decomposition and reconstruction
  const Decomposition2D<double> dec2D = haar.Wavedec(matrix, num_levels);
  const Matrix<double> matrix_rec = haar.Waverec(dec2D, matrix.Size());
  std::cout << "Error of 2D reconstruction: " << Norm(matrix.GetData() - matrix_rec.GetData()) << std::endl;

  return 0;
}
