# Simple Wavelet

Minimal C++ implementation of 1D and 2D [Wavelet transform](https://en.wikipedia.org/wiki/Wavelet_transform). The interfaces resemble [MATLAB `wavedec`](https://www.mathworks.com/help/wavelet/ref/wavedec.html) and [`wavedec2`](https://www.mathworks.com/help/wavelet/ref/wavedec2.html). This is implemented in the `Wavelet` class, which is able to compute multi-level wavelet decomposition and reconstruction of 1D and 2D data.

An auxiliary `Matrix`class is also provided to represent 2D data, whereas `std::vector` are used to represent 1D data. An additional `vector.h` header overloads standard mathematical binary and unary operators for `std::vector`.

Example usage for 1D and 2D data is provided in `main.cpp`. Run `make test` to execute the program, called `sw`. The program expects three arguments as follows
```
./sw num_rows num_cols num_levels
```
where `num_rows` and `num_cols` are two integers representing the size of a matrix, and `num_levels` the number of the wavelet decomposition levels. 

Internally the program will initialise a random matrix having size specified by the arguments and a `Wavelet` with [Haar wavlet transform filers](https://en.wikipedia.org/wiki/Haar_wavelet#Haar_transform). A column of the matrix and the full matrix will be first decomposed in `num_levels` using `Wavelet`, and then reconstructed from their decompositions. The program will print on standard output the error (l2-norm) of both reconstructions. Typical output of `make test` should be similar to
```
./sw 64 64 4
Error (l2-norm) of 1D reconstruction: 3.28881e-15
Error (l2-norm) of 2D reconstruction: 4.86914e-14
```

## Licence
MIT Licence. Copyright (c) 2018 Matteo Maggioni
