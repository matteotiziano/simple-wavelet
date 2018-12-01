#include <iostream>
#include <string>

#include "wavelet.h"

using namespace Vector;

int main(int argc, const char *argv[]) {

  // Read arguments from command line
  if(argc < 3) {
    std::cout << "usage: " << std::string(argv[0]) << " num_rows num_cols num_levels" << std::endl;
    return 1;
  }

  // TODO: handle invalid inputs

  // Size of the matrix
  const size_t num_rows = static_cast<size_t>(std::stoi(argv[1]));
  const size_t num_cols = static_cast<size_t>(std::stoi(argv[2]));
  // Number of decomposition levels
  const size_t num_levels = static_cast<size_t>(std::stoi(argv[3]));

  // Initialise a random 2D matrix using default seed
  // Note: use `matrix.Rand(std::random_device()());` to generate different
  // values at each execution
  Matrix<double> matrix(num_rows, num_cols);
  matrix.Rand();

  // Filters definition of Discrete Haar transform
  // https://en.wikipedia.org/wiki/Haar_wavelet#Haar_transform
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  const std::vector<double> Lo_D = {inv_sqrt2, inv_sqrt2};
  const std::vector<double> Hi_D = {-inv_sqrt2, inv_sqrt2};
  const std::vector<double> Lo_R = {inv_sqrt2, inv_sqrt2};
  const std::vector<double> Hi_R = {inv_sqrt2, -inv_sqrt2};

  // Instantiate a Wavelet object with Haar filters
  const Wavelet<double> haar(Lo_D, Hi_D, Lo_R, Hi_R);


  // Example of 1D decomposition and reconstruction using the first column of
  // the above matrix as input
  const std::vector<double> vector = matrix.GetColumn(0);
  const Decomposition1D<double> dec1D = haar.Wavedec(vector, num_levels);
  const std::vector<double> vector_rec = haar.Waverec(dec1D, vector.size());
  const double norm1D = Norm(vector - vector_rec);

  // Example of 2D decomposition and reconstruction
  const Decomposition2D<double> dec2D = haar.Wavedec(matrix, num_levels);
  const Matrix<double> matrix_rec = haar.Waverec(dec2D, matrix.Size());
  const double norm2D = Norm((matrix - matrix_rec).GetData());

  std::cout << "Error (l2-norm) of 1D reconstruction: " << norm1D << std::endl;
  std::cout << "Error (l2-norm) of 2D reconstruction: " << norm2D << std::endl;

  return 0;
}
