#ifndef VECTOR_H
#define VECTOR_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <vector>

namespace Vector {

  //===========================================================================
  // This header contains a definition for all the most common mathematical
  // operations on vectors implemented as overlading of standard operators.
  //===========================================================================

  //===========================================================================
  // Binary vector-vector operations
  //===========================================================================

  // Element-wise vector addition (vector1+vector2)
  template<typename T>
  std::vector<T> operator+(std::vector<T> vector1, const std::vector<T>& vector2) {
    if(vector1.size() != vector2.size()) {
      throw std::invalid_argument("std::vector operator+");
    }
    // This is actually implemented as vector1 += vector2, so the first
    // parameter is passed by copy and then overwritten with the result
    std::transform(vector1.begin(), vector1.end(), vector2.begin(), vector1.begin(), std::plus<T>());
    return vector1;
  }

  // Element-wise vector addition and assignment (vector1+=vector2)
  template<typename T>
  std::vector<T>& operator+=(std::vector<T>& vector1, const std::vector<T>& vector2) {
    if(vector1.size() != vector2.size()) {
      throw std::invalid_argument("std::vector operator+=");
    }
    // The first parameter is passed by reference, modified and returned
    std::transform(vector1.begin(), vector1.end(), vector2.begin(), vector1.begin(), std::plus<T>());
    return vector1;
  }

  // Element-wise vector subtraction (vector1-vector2)
  template<typename T>
  std::vector<T> operator-(std::vector<T> vector1, const std::vector<T>& vector2) {
    if(vector1.size() != vector2.size()) {
      throw std::invalid_argument("std::vector operator-");
    }
    std::transform(vector1.begin(), vector1.end(), vector2.begin(), vector1.begin(), std::minus<T>());
    return vector1;
  }

  // Element-wise vector subtraction and assignment (vector1-=vector2)
  template<typename T>
  std::vector<T>& operator-=(std::vector<T>& vector1, const std::vector<T>& vector2) {
    if(vector1.size()!=vector2.size()) {
      throw std::invalid_argument("std::vector operator-=");
    }
    std::transform(vector1.begin(), vector1.end(), vector2.begin(), vector1.begin(), std::minus<T>());
    return vector1;
  }

  // Element-wise vector multiplication (vector1*vector2)
  template<typename T>
  std::vector<T> operator*(std::vector<T> vector1, const std::vector<T>& vector2) {
    if(vector1.size() != vector2.size()) {
      throw std::invalid_argument("std::vector operator*");
    }
    std::transform(vector1.begin(), vector1.end(), vector2.begin(), vector1.begin(), std::multiplies<T>());
    return vector1;
  }

  // Element-wise vector multiplication and assignment (vector1*=vector2)
  template<typename T>
  std::vector<T>& operator*=(std::vector<T>& vector1, const std::vector<T>& vector2) {
    if(vector1.size() != vector2.size()) {
      throw std::invalid_argument("std::vector operator*=");
    }
    std::transform(vector1.begin(), vector1.end(), vector2.begin(), vector1.begin(), std::multiplies<T>());
    return vector1;
  }

  // Element-wise vector division (vector1/vector2)
  template<typename T>
  std::vector<T> operator/(std::vector<T> vector1, const std::vector<T>& vector2) {
    if(vector1.size() != vector2.size()) {
      throw std::invalid_argument("std::vector operator/");
    }
    std::transform(vector1.begin(), vector1.end(), vector2.begin(), vector1.begin(), std::divides<T>());
    return vector1;
  }

  // Element-wise vector division and assignment (vector1/=vector2)
  template<typename T>
  std::vector<T>& operator/=(std::vector<T>& vector1, const std::vector<T>& vector2) {
    if(vector1.size() != vector2.size()) {
      throw std::invalid_argument("std::vector operator/=");
    }
    std::transform(vector1.begin(), vector1.end(), vector2.begin(), vector1.begin(), std::divides<T>());
    return vector1;
  }

  //===========================================================================
  // Unary vector operations
  //===========================================================================

  // Element-wise vector negation (-vector)
  template<typename T>
  std::vector<T> operator-(std::vector<T> vector) {
    std::transform(vector.begin(), vector.end(), vector.begin(), std::negate<T>());
    return vector;
  }

  //===========================================================================
  // Binary vector-scalar operations
  //===========================================================================

  // Right scalar vector addition (vector+scalar)
  template<typename T>
  std::vector<T> operator+(std::vector<T> vector, T scalar) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [&scalar](T value) {
      return scalar + value;
    });
    return vector;
  }

  // Left scalar vector addition (scalar+vector)
  template<typename T>
  std::vector<T> operator+(T scalar, std::vector<T> vector) {
    return vector + scalar;
  }

  // Right scalar vector subtraction (vector-scalar)
  template<typename T>
  std::vector<T> operator-(std::vector<T> vector, T scalar) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [&scalar](T value) {
      return value - scalar;
    });
    return vector;
  }

  // Left scalar vector subtraction (scalar-vector)
  template<typename T>
  std::vector<T> operator-(T scalar, std::vector<T> vector) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [&scalar](T value) {
      return scalar - value;
    });
    return vector;
  }

  // Right scalar vector multiplication (vector*scalar)
  template<typename T>
  std::vector<T> operator*(std::vector<T> vector, T scalar) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [&scalar](T value) {
      return value * scalar;
    });
    return vector;
  }

  // Left scalar vector multiplication (scalar*vector)
  template<typename T>
  std::vector<T> operator*(T scalar, std::vector<T> vector) {
    return vector * scalar;
  }

  // Right scalar vector division (vector/scalar)
  template<typename T>
  std::vector<T> operator/(std::vector<T> vector, T scalar) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [&scalar](T value) {
      return value / scalar;
    });
    return vector;
  }

  // Left scalar vector division (scalar/vector)
  template<typename T>
  std::vector<T> operator/(T scalar, std::vector<T> vector) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [&scalar](T value) {
      return scalar / value;
    });
    return vector;
  }

  // Vector power to the given scalar (vector^scalar)
  template<typename T>
  std::vector<T> operator^(std::vector<T> vector, T scalar) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [&scalar](T value) {
      return std::pow(value, scalar);
    });
    return vector;
  }

  //===========================================================================
  // Utils
  //===========================================================================

  // Element-wise round of vector
  template<typename T>
  void Round(std::vector<T>& vector) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [](T value) {
      return std::round(value);
    });
  }

  // Element-wise floor of vector
  template<typename T>
  void Floor(std::vector<T>& vector) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [](T value) {
      return std::floor(value);
    });
  }

  // Element-wise ceil of vector
  template<typename T>
  void Ceil(std::vector<T>& vector) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [](T value) {
      return std::ceil(value);
    });
  }

  // Element-wise clamp of vector
  template<typename T>
  void Clamp(std::vector<T>& vector, T min_value, T max_value) {
    std::transform(vector.begin(), vector.end(), vector.begin(), [&min_value, &max_value](T value) {
      return std::min(std::max(value, min_value), max_value);
    });
  }

  // Get maximum value of vector
  template<typename T>
  T Max(const std::vector<T>& vector) {
    return *std::max_element(vector.begin(), vector.end());
  }

  // Get minimum value of vector
  template<typename T>
  T Min(const std::vector<T>& vector) {
    return *std::min_element(vector.begin(), vector.end());
  }

  // Sum of elements in vector
  template<typename T>
  T Sum(const std::vector<T>& vector) {
    T val = 0;
    val = std::accumulate(vector.begin(), vector.end(), val);
    return val;
  }

  // Product of elements in vector
  template<typename T>
  T Prod(const std::vector<T>& vector) {
    T val = 1;
    val = std::accumulate(vector.begin(), vector.end(), val, std::multiplies<T>());
    return val;
  }

  // Dot-product between vectors
  template<typename T>
  T Dot(const std::vector<T>& vector1, const std::vector<T>& vector2) {
    return Sum(vector1 * vector2);
  }

  // Lp-norm of vector
  template<typename T>
  T LpNorm(const std::vector<T>& vector, T p) {
    return Sum(vector ^ p) ^ (static_cast<T>(1) / p);
  }

  // L2-norm of vector
  template<typename T>
  T Norm(const std::vector<T>& vector) {
    return std::sqrt(Dot(vector, vector));
  }

  // Normalise vector
  template<typename T>
  std::vector<T> Normalise(const std::vector<T>& vector) {
    return vector / Norm(vector);
  }

  // Convolution types
  enum ConvolutionType {
    kFullConvolution = 0,
    kSameConvolution,
    kValidConvolution
  };

  // Convolution of vector with given filter
  template<typename T>
  std::vector<T> Conv(const std::vector<T>& vector, const std::vector<T>& filter, ConvolutionType convolutionType = kFullConvolution) {
    const size_t vector_size = vector.size();
    const size_t filter_size = filter.size();
    if (filter_size > vector_size) {
      throw std::invalid_argument("Conv");
    }
    // Compute convolution
    size_t minI = 0;
    size_t maxI = vector_size + filter_size - 1;
    if (convolutionType == kSameConvolution) {
      minI = filter_size / 2;
      maxI = minI + vector_size;
    }
    else if (convolutionType == kValidConvolution) {
      minI = filter_size - 1;
      maxI = vector_size;
    }
    std::vector<T> output;
    output.reserve(maxI - minI);
    for (size_t i = minI; i < maxI; i++) {
      const size_t minJ = (i < filter_size - 1) ? 0 : i - filter_size + 1;
      const size_t maxJ = (i >= vector_size - 1) ? vector_size : i + 1;
      T value = static_cast<T>(0);
      for (size_t j = minJ; j < maxJ; j++) {
        value += vector[j] * filter[i - j];
      }
      output.push_back(value);
    }
    return output;
  }

} // namespace Vector

//=============================================================================
// Output to std::ostream
//=============================================================================

// Overloaded operator<< for std::vector<T>
template<typename T>
inline std::ostream& operator<<(std::ostream& out, const std::vector<T>& vector) {
  for(size_t i = 0, endI = vector.size(); i < endI; i++) {
    out << std::setw(8) << std::setprecision(4) << std::left << vector.at(i);
  }
  return out;
}

#endif // VECTOR_H
