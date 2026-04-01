HYLU----Hybrid Parallel Sparse LU Factorization
=========
HYLU is a general-purpose parallel solver designed for efficiently solving sparse linear systems ($\bf{Ax}=\bf b$) on multi-core shared-memory machines. It employs an innovative parallel up-looking matrix factorization algorithm, which dynamically adapts to varying matrix sparsity patterns by leveraging hybrid numerical kernels.


HYLU delivers high-performance matrix factorization for large-scale sparse linear systems from multiple engineering and scientific domains, including circuit simulation, power systems, computational fluid dynamics (CFD), electromagnetics, and structural analysis. The solver efficiently handles linear systems arising from finite element analysis, 2D/3D modeling, and optimization problems. Unsymmetric, symmetric indefinite, and symmetric positive-definite matrices are supported.


Performance Results
============

Please see [doc/results.pdf](https://github.com/chenxm1986/hylu/blob/main/doc/results.pdf) for details. 

For one-time solve, HYLU is 1.70X faster than Intel MKL PARDISO on geometric mean in the total time of preprocessing, factorization, and solve.

For repeated solve, HYLU is 2.53X faster than Intel MKL PARDISO on geometric mean in the time of factorization and solve.

For the solution accuracy, HYLU achieves an order of magnitude higher accuracy than Intel MKL PARDISO.


System Requirements
=============
+ Hardware requirement: X64 CPU with AVX2 and FMA instructions supported

+ Software requirement: Windows 7/10/11 or CentOS 7/8 (mainstream Linux distributions are also supported)


Notes on Library and Integer Bitwidths
============
Only x64 libraries are provided. This means that, a 64-bit Windows or Linux operating system is needed.

Functions for both 32-bit integers and 64-bit integers are provided. The latter has '_L' in the function names. The integer bitwidth only limits the size of the input matrix. The internal data structures always use 64-bit integers.


History
============
+ Version 20260331
	+ Added support for symmetric indefinite and symmetric positive-definite matrices
	+ Removed requirement of ensuring non-zero diagonals when using custom ordering

+ Version 20260205
	+ Added functions to calculate determinant and condition number
	+ Added solve of conjugate transposed system (only for complex number system)

+ Version 20260126
	+ Updated parm[23] (whether to ensure consistent symbolic results between sequential and parallel executions)

+ Version 20260116
	+ Performance improvement

+ Version 20260104
	+ Performance improvement for factorization

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