HYLU----Hybrid Parallel Sparse LU Factorization
=========
HYLU is a general-purpose parallel solver designed for efficiently solving sparse linear systems ($\bf{Ax}=\bf b$) on multi-core shared-memory machines. It employs an innovative parallel up-looking LU factorization algorithm, which dynamically adapts to varying matrix sparsity patterns by leveraging hybrid numerical kernels.


HYLU delivers high-performance LU factorization for large-scale sparse linear systems from multiple engineering and scientific domains, including circuit simulation, power systems, computational fluid dynamics (CFD), electromagnetics, and structural analysis. The solver efficiently handles linear systems arising from finite element analysis, 2D/3D modeling, and optimization problems.


Performance Results
============
For a wide range of matrices with different sparsities, HYLU achieves a **2.04X** speedup on geometric mean in numerical factorization compared with Intel MKL PARDISO, while the preprocessing and forward-backward substitution phases are also faster (**1.56X** and **1.46X** speedups, respectively). HYLU offers an optimization option for repeated solve of linear systems with an identical sparse pattern in the coefficient matrix. In this scenario, HYLU achieves a **2.58X** geometric mean speedup in numerical factorization over Intel MKL PARDISO, while the forward-backward substitution phase is slightly faster (**1.35X** speedup). Please see [doc/results.pdf](https://github.com/chenxm1986/hylu/blob/main/doc/results.pdf) for details. Generally, HYLU is much faster than Intel MKL PARDISO for highly sparse matrices, while for relatively dense matrices, HYLU achieves similar performance to Intel MKL PARDISO.




Notes on Library and Integer Bitwidths
============
Only x64 libraries are provided. This means that, a 64-bit Windows or Linux operating system is needed.

Functions for both 32-bit integers and 64-bit integers are provided. The latter has '_L' in the function names. The integer bitwidth only limits the size of the input matrix. The internal data structures always use 64-bit integers.


History
============
+ Version 20251222
	+ Updated parallel nested dissection ordering

+ Version 20251212
	+ Added a function to solve multiple right-hand-side vectors
	+ Updated nested dissection for repeated solving

+ Version 20251203
	+ Added a function to support user-provided ordering
	+ Added complex number support

+ Version 20251015
	+ Initial release

+ Version before 20251015
	+ Test version


Publications
============
[1] Xiaoming Chen, "HYLU: Hybrid Parallel Sparse LU Factorization", arXiv: 2509.07690, https://arxiv.org/abs/2509.07690.

Author
============
Please visit [Xiaoming Chen's personal page](http://people.ucas.edu.cn/~chenxm).