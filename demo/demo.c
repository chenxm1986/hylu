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
#pragma comment(lib, "hylu.lib")
#endif

static bool ReadMtxFile(const char file[], int *n, int **ap, int **ai, double **ax)
{
    FILE *fp = fopen(file, "r");
    if (NULL == fp)
    {
        printf("Cannot open file \"%s\".\n", file);
        return false;
    }

    char buf[256] = "\0";
    bool first = true;
    int pc = 0;
    int ptr = 0;
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
                int r, c, nz;
                sscanf(p, "%d %d %d", &r, &c, &nz);
                if (r != c)
                {
                    printf("Matrix is not square because row = %d and column = %d.\n", r, c);
                    fclose(fp);
                    return false;
                }

                *n = r;
                *ap = (int *)malloc(sizeof(int) * (*n + 1));
                *ai = (int *)malloc(sizeof(int) * nz);
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
                int r, c;
                double v;
                sscanf(p, "%d %d %lf", &r, &c, &v);
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

static double L1NormOfResidual(const int n, const int ap[], const int ai[], const double ax[], const double x[], const double b[], bool row0_col1)
{
    if (row0_col1)
    {
        double *bb = (double *)malloc(sizeof(double) * n);
        memcpy(bb, b, sizeof(double) * n);
        for (int i = 0; i < n; ++i)
        {
            const double xx = x[i];
            const int start = ap[i];
            const int end = ap[i + 1];
            for (int p = start; p < end; ++p)
            {
                bb[ai[p]] -= xx * ax[p];
            }
        }
        double s = 0., s2 = 0.;
        for (int i = 0; i < n; ++i)
        {
            s += fabs(bb[i]);
            s2 += fabs(b[i]);
        }
        free(bb);
        return 0. == s2 ? 0. : s / s2;
    }
    else
    {
        double s = 0., s2 = 0.;
        for (int i = 0; i < n; ++i)
        {
            double r = b[i];
            const int start = ap[i];
            const int end = ap[i + 1];
            for (int p = start; p < end; ++p)
            {
                const int j = ai[p];
                r -= ax[p] * x[j];
            }
            s += fabs(r);
            s2 += fabs(b[i]);
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
        printf("Usage: ./demo <mtx file> <# of threads>\n");
        printf("Example: ./demo ss1.mtx 4\n");
        printf("mtx files can be downloaded from https://sparse.tamu.edu\n");
        return -1;
    }

    int n;
    int *ap = NULL;
    int *ai = NULL;
    double *ax = NULL;
    double *b = NULL, *x = NULL;
    void *instance = NULL;
    long long *parm = NULL;
    int ret = 0;
    double res;
    unsigned int rd = 1234;
    double mantissa, exponent;
    double cond;

    if (!ReadMtxFile(argv[1], &n, &ap, &ai, &ax)) goto RETURN;
    printf("N(A) = %d.\n# NZ(A) = %d.\n", n, ap[n]);

    b = (double *)malloc(sizeof(double) * n * 2);
    if (NULL == b)
    {
        printf("Malloc for b and x failed.\n");
        goto RETURN;
    }
    x = b + n;
    for (int i = 0; i < n; ++i)
    {
        b[i] = (double)Rand(&rd) / (double)0xFFFFFFFF * 100. - 50.;
        x[i] = 0.;
    }

    ret = HYLU_CreateSolver(&instance, &parm, atoi(argv[2]));
    if (ret < 0)
    {
        printf("Failed to create solver instance, return code = %d.\n", ret);
        goto RETURN;
    }
    parm[1] = 1;

    ret = HYLU_Analyze(instance, false, n, ap, ai, ax);
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
    printf("# Supernodes = %d.\n", (int)parm[9]);

    ret = HYLU_Factorize(instance, ax);
    if (ret < 0)
    {
        printf("Numerical factorization failed, return code = %d.\n", ret);
        goto RETURN;
    }
    printf("Numerical factorization time = %g.\n", parm[7] * 1.e-6);
    printf("# Swapped pivots = %d.\n", (int)parm[8]);
    printf("# Perturbed pivots = %d.\n", (int)parm[11]);

    ret = HYLU_Solve(instance, false, b, x);
    if (ret < 0)
    {
        printf("Triangular solve failed, return code = %d.\n", ret);
        goto RETURN;
    }
    printf("Triangular solve time = %g.\n", parm[7] * 1.e-6);
    printf("# Refinements = %d.\n", (int)parm[16]);

    res = L1NormOfResidual(n, ap, ai, ax, x, b, false);
    printf("Residual (||Ax-b||/||b||) = %g.\n", res);

    ret = HYLU_Determinant(instance, &mantissa, &exponent);
    if (ret < 0)
    {
        printf("Determinant calculation failed, return code = %d.\n", ret);
        goto RETURN;
    }
    printf("Determinant = %g * 10**(%.0lf).\n", mantissa, exponent);

    ret = HYLU_ConditionNumber(instance, &cond);
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
    free(b);
    HYLU_DestroySolver(instance);
	return 0;
}