/*
* HYLU (Hybrid Parallel Sparse LU Factorization) is a general-purpose parallel solver designed for efficiently solving sparse linear systems (Ax=b) 
* on multi-core shared-memory machines.
*/

#ifndef __HYLU_H__
#define __HYLU_H__

/********** error code **********
*  0:   successful
* -1:   invalid instance handle
* -2:   argument error
* -3:   invalid matrix
* -4:   out of memory
* -5:   structurally singular
* -6:   numerically singular
* -7:   threads error
* -8:   calling procedure error
* -9:   integer overflow
* -10:  internal error
********************************/

/********** parm spec **********
* parm[0]: output, version
* parm[1]: input, timer. [default 0]: no timer | >0: microsecond-level timer | <0: millisecond-level timer
* parm[2]: input, ordering method. [default 0]: automatic select | 1: approximate minimum degree | 2: approximate minimum degree variant
		   | 3: nested dissection | 4: nested dissection variant | 5: best of 1 and 2 | 6: best of 3 and 4 | 7: best of 1-4
* parm[3]: input, ordering method switch point. [default 0]: automatic control
* parm[4]: output, selected ordering method (1-4)
* parm[5]: input, minimum # of supernode columns. [default 32]
* parm[6]: input, maximum # of supernode rows. [default 0]: automatic control
* parm[7]: output, time (in microsecond) of last function call
* parm[8]: output, # of off-diagonal pivots
* parm[9]: output, # of supernodes
* parm[10]: input, pivot perturbation, zero or small pivot will be replaced by (10**parm[10])*||A||. [default -15]
* parm[11]: output, # of perturbed pivots
* parm[12]: output, current memory usage (in bytes), or needed memory size (in bytes) when -4 is returned
* parm[13]: output, maximum memory usage (in bytes)
* parm[14]: output, thread # stored in 3 shorts, from lowest to highest: physical core # (may be incorrect), logical core #, created threads #. highest short is not used
* parm[15]: input, maximum # of refinement iterations. [default 0]: automatic control | >0: perform at most parm[15] iterations | <0: force to perform -parm[15] iterations
* parm[16]: output, # of refinement iterations performed
* parm[17]: output, # of nonzeros in L, including diagonal
* parm[18]: output, # of nonzeros in U, excluding diagonal
* parm[19]: output, # of flops of factorization (excluding scaling)
* parm[20]: output, # of flops of solving (excluding scaling)
* parm[21]: input, whether to scale matrix. [default >0]: dynamic scaling | <0: static scaling | 0: no scaling
* parm[22]: input, symbolic factorization method. [default 0]: automatic control | >0: unsymmetric symbolic factorization | <0: symmetric symbolic factorization
* parm[23]: input, whether to use parallel nested dissection. [default 1]: enabled | 0: disabled
********************************/

#ifndef __cplusplus
/*
* stdbool.h is added in C99. If the file does not exist, simply replace the following line with "#define bool char"
*/
#include <stdbool.h>
#endif

#define _IN_
#define _OUT_

#ifdef __cplusplus
extern "C"
{
#endif

/*
* Creates solver instance, retrieves parameter array, and creates threads
* @instance: will return pointer to solver instance
* @parm: will return pointer to parameter array; do not free this array; set to NULL if not needed
* @threads: # of threads to be created; cannot exceed # of logical cores
*/
int HYLU_CreateSolver
(
	_OUT_ void **instance,
	_OUT_ long long **parm, /*can be NULL if not needed*/
	_IN_ int threads /*1: sequential. 0: use all physical cores. -1: use all logical cores*/
);
int HYLU_L_CreateSolver
(
	_OUT_ void **instance,
	_OUT_ long long **parm, /*can be NULL if not needed*/
	_IN_ int threads /*1: sequential. 0: use all physical cores. -1: use all logical cores*/
);

/*
* Frees memory, terminates threads, and destroys solver instance
* @instance: solver instance
*/
int HYLU_DestroySolver
(
	_IN_ void *instance
);
int HYLU_L_DestroySolver
(
	_IN_ void *instance
);

/*
* Analyzes matrix for ordering and symbolic factorization
* @instance: solver instance
* @repeat: whether linear systems will be solved repeatedly with identical matrix structure (false for one-time solving)
* @n: matrix dimension
* @ap: integer array of length n+1, matrix row pointers
* @ai: integer array of length ap[n], matrix column indexes
* @ax: double array of length ap[n], matrix values
*/
int HYLU_Analyze
(
	_IN_ void *instance,
	_IN_ bool repeat,
	_IN_ int n,
	_IN_ const int ap[],
	_IN_ const int ai[],
	_IN_ const double ax[] /*can be NULL (static pivoting and static scaling will be disabled)*/
);
int HYLU_L_Analyze
(
	_IN_ void *instance,
	_IN_ bool repeat,
	_IN_ long long n,
	_IN_ const long long ap[],
	_IN_ const long long ai[],
	_IN_ const double ax[] /*can be NULL (static pivoting and static scaling will be disabled)*/
);

/*
* Factorizes matrix with diagonal block pivoting
* @instance: solver instance
* @ax: double array of length ap[n], matrix values
*/
int HYLU_Factorize
(
	_IN_ void *instance,
	_IN_ const double ax[] /*can be different from that used for HYLU_Analyze*/
);
int HYLU_L_Factorize
(
	_IN_ void *instance,
	_IN_ const double ax[] /*can be different from that used for HYLU_L_Analyze*/
);

/*
* Solves Ax=b after A is factorized
* @instance: solver instance
* @transpose: whether to solve (A**T)x=b, note that HYLU uses row-major order by default
* @b: double array of length n to specify right-hand-side vector
* @x: double array of length n to get solution
*/
int HYLU_Solve
(
	_IN_ void *instance,
	_IN_ bool transpose, /*false for row mode, true for column mode*/
	_IN_ const double b[],
	_OUT_ double x[] /*x space can overlap b space*/
);
int HYLU_L_Solve
(
	_IN_ void *instance,
	_IN_ bool transpose, /*false for row mode, true for column mode*/
	_IN_ const double b[],
	_OUT_ double x[] /*x space can overlap b space*/
);

#ifdef __cplusplus
}
#endif

#endif