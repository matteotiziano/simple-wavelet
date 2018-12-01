# Simple Wavelet

Minimal C++ implementation of 1D and 2D [Wavelet transform](https://en.wikipedia.org/wiki/Wavelet_transform). The interfaces resemble [MATLAB `wavedec`](https://www.mathworks.com/help/wavelet/ref/wavedec.html) and [`wavedec2`](https://www.mathworks.com/help/wavelet/ref/wavedec2.html). This is implemented in the `Wavelet` class, which is able to compute multi-level wavelet decomposition and reconstruction of 1D and 2D data.

An auxiliary `Matrix`class is also provided to represent 2D data, whereas `std::vector` are used to represent 1D data. An additional `vector.h` header overloads standard mathematical binary and unary operators for `std::vector`.

Example usage for 1D and 2D data is provided in `main.cpp`. Run `make` to compile and execute the program, called `sw`. The program expects three arguments as follows
```
./sw num_rows num_cols num_levels
```
which will be used to initialise a random matrix having size `num_rows`x`num_cols`, and a `Wavelet` transform with [Haar transform filers](https://en.wikipedia.org/wiki/Haar_wavelet#Haar_transform) to compute a 1D and 2D decomposition with `num_levels`.

The matrix will be used as input for the 2D decomposition, whereas the first column of the matrix (having length `num_rows`) will be used as 1D data. Both 1D and 2D data will be reconstructed from their decomposition, and the reconstruction error (l2-norm) will be finally printed to standard output. Note that the seed for the random value generation is fixed, thus the results will be
```
./sw 64 64 4
Error (l2-norm) of 1D reconstruction: 3.30226e-15
Error (l2-norm) of 2D reconstruction: 5.0788e-14
```

## Licence
MIT Licence. Copyright (c) 2018 Matteo Maggioni
