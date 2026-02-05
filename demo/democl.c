/*this demo code reads a real number matrix and generates the imaginary parts randomly*/
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS		1
#define _CRT_NONSTDC_NO_WARNINGS	1
#define _CRT_SECURE_NO_DEPRECATE	1
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "hylu.h"
#ifdef _MSC_VER
#pragma comment(lib, "hylu_cl.lib")
#endif

static bool ReadMtxFile(const char file[], long long *n, long long **ap, long long **ai, double **ax)
{
    FILE *fp = fopen(file, "r");
    if (NULL == fp)
    {
        printf("Cannot open file \"%s\".\n", file);
        return false;
    }

    char buf[256] = "\0";
    bool first = true;
    long long pc = 0;
    long long ptr = 0;
    while (fgets(buf, 256, fp) != NULL)
    {
        const char *p = buf;
        while (*p != '\0')
        {
            if (' ' == *p || '\t' == *p || '\r' == *p || '\n' == *p) ++p;
            else break;
        }

        if (*p == '\0') continue;
        else if (*p == '%') continue;
        else
        {
            if (first)
            {
                first = false;
                long long r, c, nz;
                sscanf(p, "%lld %lld %lld", &r, &c, &nz);
                if (r != c)
                {
                    printf("Matrix is not square because row = %lld and column = %lld.\n", r, c);
                    fclose(fp);
                    return false;
                }

                *n = r;
                *ap = (long long *)malloc(sizeof(long long) * (*n + 1));
                *ai = (long long *)malloc(sizeof(long long) * nz);
                *ax = (double *)malloc(sizeof(double) * nz);
                if (NULL == *ap || NULL == *ai || NULL == *ax)
                {
                    printf("Malloc for matrix failed.\n");
                    fclose(fp);
                    return false;
                }
                (*ap)[0] = 0;
            }
            else
            {
                long long r, c;
                double v;
                sscanf(p, "%lld %lld %lf", &r, &c, &v);
                --r;
                --c;
                (*ai)[ptr] = r;
                (*ax)[ptr] = v;
                if (c != pc)
                {
                    (*ap)[c] = ptr;
                    pc = c;
                }
                ++ptr;
            }
        }
    }
    (*ap)[*n] = ptr;

    fclose(fp);
    return true;
}

static __inline void __cmulsub(complex_t z, const complex_t a, const complex_t b) /* z-=a*b */
{
	const double a0 = a[0];
	const double a1 = a[1];
	const double b0 = b[0];
	const double b1 = b[1];
	z[0] -= a0 * b0 - a1 * b1;
	z[1] -= a0 * b1 + a1 * b0;
}

static __inline void __cmulsub_conj(complex_t z, const complex_t a, const complex_t b) /* z-=conj(a)*b */
{
	const double a0 = a[0];
	const double a1 = a[1];
	const double b0 = b[0];
	const double b1 = b[1];
	z[0] -= a0 * b0 + a1 * b1;
	z[1] -= a0 * b1 - a1 * b0;
}

static double L1NormOfResidual(const long long n, const long long ap[], const long long ai[], const complex_t ax[], const complex_t x[], const complex_t b[], char mode)
{
    if (mode)
    {
        complex_t *bb = (complex_t *)malloc(sizeof(complex_t) * n);
        memcpy(bb, b, sizeof(complex_t) * n);
        if (mode > 0)
        {
            for (long long i = 0; i < n; ++i)
            {
                const complex_t xx = { x[i][0], x[i][1] };
                const long long start = ap[i];
                const long long end = ap[i + 1];
                for (long long p = start; p < end; ++p)
                {
                    __cmulsub(bb[ai[p]], ax[p], xx);
                }
            }
        }
        else
        {
            for (long long i = 0; i < n; ++i)
            {
                const complex_t xx = { x[i][0], x[i][1] };
                const long long start = ap[i];
                const long long end = ap[i + 1];
                for (long long p = start; p < end; ++p)
                {
                    __cmulsub_conj(bb[ai[p]], ax[p], xx);
                }
            }
        }
        double s = 0., s2 = 0.;
        for (long long i = 0; i < n; ++i)
        {
            s += hypot(bb[i][0], bb[i][1]);
            s2 += hypot(b[i][0], b[i][1]);
        }
        free(bb);
        return 0. == s2 ? 0. : s / s2;
    }
    else
    {
        double s = 0., s2 = 0.;
        for (long long i = 0; i < n; ++i)
        {
            complex_t r = { b[i][0], b[i][1]};
            const long long start = ap[i];
            const long long end = ap[i + 1];
            for (long long p = start; p < end; ++p)
            {
                const long long j = ai[p];
                __cmulsub(r, ax[p], x[j]);
            }
            s += hypot(r[0], r[1]);
            s2 += hypot(b[i][0], b[i][1]);
        }
        return 0. == s2 ? 0. : s / s2;
    }
}

static __inline unsigned int Rand(unsigned int *x)
{
    *x = *x * 134775813 + 1;
    return *x;
}

int main(int argc, const char *argv[])
{
    if (argc < 3)
    {
        printf("Usage: ./democl <mtx file> <# of threads>\n");
        printf("Example: ./democl ss1.mtx 4\n");
        printf("mtx files can be downloaded from https://sparse.tamu.edu\n");
        printf("This demo code reads a real number matrix and generates the imaginary parts randomly\n");
        return -1;
    }

    long long n, nz;
    long long *ap = NULL;
    long long *ai = NULL;
    double *ax = NULL;
    complex_t *cx = NULL;
    complex_t *b = NULL, *x = NULL;
    void *instance = NULL;
    long long *parm = NULL;
    int ret = 0;
    double res;
    unsigned int rd = 1234;
    complex_t mantissa;
    double exponent;
    double cond;

    if (!ReadMtxFile(argv[1], &n, &ap, &ai, &ax)) goto RETURN;
    printf("N(A) = %lld.\n# NZ(A) = %lld.\n", n, ap[n]);

    nz = ap[n];
    cx = (complex_t *)malloc(sizeof(complex_t) * nz);
    if (NULL == cx)
    {
        printf("Malloc for cx failed.\n");
        goto RETURN;
    }
    for (long long i = 0; i<nz; ++i)
    {
        cx[i][0] = ax[i];
        cx[i][1] = ax[i] * ((double)Rand(&rd) / (double)0xFFFFFFFF * 8. - 4.);
    }

    b = (complex_t *)malloc(sizeof(complex_t) * n * 2);
    if (NULL == b)
    {
        printf("Malloc for b and x failed.\n");
        goto RETURN;
    }
    x = b + n;
    for (long long i = 0; i < n; ++i)
    {
        b[i][0] = (double)Rand(&rd) / (double)0xFFFFFFFF * 100. - 50.;
        b[i][1] = (double)Rand(&rd) / (double)0xFFFFFFFF * 100. - 50.;
        x[i][0] = 0.;
        x[i][1] = 0.;
    }

    ret = HYLU_CL_CreateSolver(&instance, &parm, atoi(argv[2]));
    if (ret < 0)
    {
        printf("Failed to create solver instance, return code = %d.\n", ret);
        goto RETURN;
    }
    parm[1] = 1;

    ret = HYLU_CL_Analyze(instance, false, n, ap, ai, cx);
    if (ret < 0)
    {
        printf("Matrix analysis failed, return code = %d.\n", ret);
        goto RETURN;
    }
    printf("Matrix analysis time = %g.\n", parm[7] * 1.e-6);
    printf("# NZ(L) = %lld.\n", parm[17]);
    printf("# NZ(U) = %lld.\n", parm[18]);
    printf("# FLOPS(factor) = %lld.\n", parm[19]);
    printf("# FLOPS(solve) = %lld.\n", parm[20]);
    printf("# Supernodes = %lld.\n", parm[9]);

    ret = HYLU_CL_Factorize(instance, cx);
    if (ret < 0)
    {
        printf("Numerical factorization failed, return code = %d.\n", ret);
        goto RETURN;
    }
    printf("Numerical factorization time = %g.\n", parm[7] * 1.e-6);
    printf("# Swapped pivots = %lld.\n", parm[8]);
    printf("# Perturbed pivots = %lld.\n", parm[11]);

    ret = HYLU_CL_Solve(instance, 0, b, x);
    if (ret < 0)
    {
        printf("Triangular solve failed, return code = %d.\n", ret);
        goto RETURN;
    }
    printf("Triangular solve time = %g.\n", parm[7] * 1.e-6);
    printf("# Refinements = %d.\n", (int)parm[16]);

    res = L1NormOfResidual(n, ap, ai, cx, x, b, 0);
    printf("Residual (||Ax-b||/||b||) = %g.\n", res);

    ret = HYLU_CL_Determinant(instance, &mantissa, &exponent);
    if (ret < 0)
    {
        printf("Determinant calculation failed, return code = %d.\n", ret);
        goto RETURN;
    }
    printf("Determinant = (%g, %g) * 10**(%.0lf).\n", mantissa[0], mantissa[1], exponent);

    ret = HYLU_CL_ConditionNumber(instance, &cond);
    if (ret < 0)
    {
        printf("Condition number calculation failed, return code = %d.\n", ret);
        goto RETURN;
    }
    printf("Condition number = %g.\n", cond);

RETURN:
    free(ap);
    free(ai);
    free(ax);
    free(cx);
    free(b);
    HYLU_CL_DestroySolver(instance);
	return 0;
}