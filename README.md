HYLU----Hybrid Parallel Sparse LU Factorization
=========
HYLU is a general-purpose parallel solver designed for efficiently solving sparse linear systems ($\bf{Ax}=\bf b$) on multi-core shared-memory machines. It employs an innovative parallel up-looking LU factorization algorithm, which dynamically adapts to varying matrix sparsity patterns by leveraging hybrid numerical kernels.

Performance Results
============
For a wide range of matrices with different sparsities, HYLU achieves a 1.74X speedup on geometric mean in numerical factorization compared with Intel MKL PARDISO, while the preprocessing and forward-backward substitution performance are both similar (1.05X and 1.22X, respectively). Please see [doc/results.pdf](https://github.com/chenxm1986/hylu/blob/main/doc/results.pdf) for details.

Notes on Library and Integer Bitwidths
============
Only x64 libraries are provided. This means that, a 64-bit Windows or Linux operating system is needed.

Functions for both 32-bit integers and 64-bit integers are provided. The latter has '_L' in the function names. The integer bitwidth only limits the size of the input matrix. The internal data structures always use 64-bit integers.

Author
============
Please visit [Xiaoming Chen's personal page](http://people.ucas.edu.cn/~chenxm).