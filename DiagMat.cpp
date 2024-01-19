#include <iostream>

#include "DiagMat.h"

// external C library.
extern "C"
{
    // Declaration of zgeev from LAPACK
    extern void zgeev_(char *jobvl, char *jobvr, int *n, std::complex<double> *a, int *lda, std::complex<double> *w,
                       std::complex<double> *vl, int *ldvl, std::complex<double> *vr, int *ldvr, std::complex<double> *work,
                       int *lwork, double *rwork, int *info);
}
//
//
void diagGenComplexMat(int &dim,
                       std::complex<double> **mat,
                       std::complex<double> **pmat,
                       std::complex<double> *eigenvalue,
                       std::complex<double> **pmatinv,
                       int &info)
{

    std::complex<double> *a = 0;
    a = new std::complex<double>[dim * dim];
    //
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            a[i * dim + j] = mat[i][j];
        }
    }

    char jobvl = 'N';
    char jobvr = 'V';
    int n = dim;
    int lda = dim;
    int ldvl = dim;
    int ldvr = dim;
    int lwork_query = -1; // Query optimal workspace size
    // int info;

    // Variables for workspace query
    std::complex<double> work_query;
    double rwork_query;

    // Call zgeev with query to get optimal workspace size
    // zgeev_()
    zgeev_(&jobvl, &jobvr, &n, a, &lda, nullptr, nullptr, &ldvl, nullptr, &ldvr, &work_query, &lwork_query, &rwork_query, &info);

    // Check for successful query
    if (info != 0)
    {
        std::cerr << "Error: zgeev workspace query failed with info = " << info << std::endl;
        exit(1);
    }

    // Get optimal workspace size
    int lwork = static_cast<int>(work_query.real()) + 1; // Add 1 as a safety margin

    // Allocate arrays for zgeev
    std::complex<double> *w = new std::complex<double>[n];
    std::complex<double> *vr = new std::complex<double>[ldvr * n];
    std::complex<double> *work = new std::complex<double>[lwork];
    double *rwork = new double[2 * n];

    // Call zgeev with allocated workspace
    zgeev_(&jobvl, &jobvr, &n, a, &lda, w, nullptr, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

    // Check for successful computation
    if (info != 0)
    {
        std::cerr << "Error: zgeev failed with info = " << info << std::endl;
        exit(1);
    }

    // Print eigenvalues
    std::cout << "Eigenvalues:" << std::endl;
    for (int i = 0; i < dim; ++i)
    {
        std::cout << w[i] << std::endl;
    }
}
