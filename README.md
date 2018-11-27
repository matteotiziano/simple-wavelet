# Simple Wavelet

Minimal C++ implementation of 1D and 2D [Wavelet transform](https://en.wikipedia.org/wiki/Wavelet_transform). The interfaces resemble [MATLAB `wavedec`](https://www.mathworks.com/help/wavelet/ref/wavedec.html) and [`wavedec2`](https://www.mathworks.com/help/wavelet/ref/wavedec2.html).

An auxiliary `Matrix`class is also provided to represent 2D data, whereas `std::vector` are used to represent 1D data. An additional `vector.h` header overloads standard mathematical binary and unary operators for `std::vector`.

Example usage is provided in `main.cpp`. Run `make test` to execute the program.

## Licence
MIT Licence. Copyright (c) 2018 Matteo Maggioni
