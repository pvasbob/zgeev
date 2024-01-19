#include <iostream>
#include <complex>
#include <vector>
// #include <lapacke.h>

extern "C"
{
    // Declaration of zgeev from LAPACK
    extern void zgeev_(char *jobvl, char *jobvr, int *n, std::complex<double> *a, int *lda, std::complex<double> *w,
                       std::complex<double> *vl, int *ldvl, std::complex<double> *vr, int *ldvr, std::complex<double> *work,
                       int *lwork, double *rwork, int *info);
}

int main()
{
    // Constants
    const int dim = 3;

    // Matrix definition
    // std::vector<std::complex<double>> mat = {
    // {0.0, 0.0}, {1.0, -1.0}, {0.0, 0.0}, {1.0, 1.0}, {0.0, 0.0}, {1.0, -1.0}, {0.0, 0.0}, {1.0, 1.0}, {0.0, 0.0}};
    //

    std::complex<double> **mat;
    //
    mat = new std::complex<double> *[3];
    for (int i = 0; i < 3; i++)
    {
        mat[i] = new std::complex<double>[3];
    };

    mat[0][0] = std::complex<double>(0.0, 0.0);
    mat[0][1] = std::complex<double>(1.0, -1.0);
    mat[0][2] = std::complex<double>(0.0, 0.0);
    //
    mat[1][0] = std::complex<double>(1.0, 1.0);
    mat[1][1] = std::complex<double>(0.0, 0.0);
    mat[1][2] = std::complex<double>(1.0, -1.0);
    //
    mat[2][0] = std::complex<double>(0.0, 0.0);
    mat[2][1] = std::complex<double>(1.0, 1.0);
    mat[2][2] = std::complex<double>(0.0, 0.0);

    // mat[0][0] = std::complex<double>(0.0, 0.0);
    // mat[0][1] = std::complex<double>(1.0, 1.0);
    // mat[0][2] = std::complex<double>(0.0, 0.0);
    ////
    // mat[1][0] = std::complex<double>(1.0, -1.0);
    // mat[1][1] = std::complex<double>(0.0, 0.0);
    // mat[1][2] = std::complex<double>(1.0, 1.0);
    ////
    // mat[2][0] = std::complex<double>(0.0, 0.0);
    // mat[2][1] = std::complex<double>(1.0, -1.0);
    // mat[2][2] = std::complex<double>(0.0, 0.0);

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

    //
    // a[0] = std::complex<double>(0.0, 0.0);
    // a[1] = std::complex<double>(1.0, -1.0);
    // a[2] = std::complex<double>(0.0, 0.0);
    ////
    // a[3] = std::complex<double>(1.0, 1.0);
    // a[4] = std::complex<double>(0.0, 0.0);
    // a[5] = std::complex<double>(1.0, -1.0);
    ////
    // a[6] = std::complex<double>(0.0, 0.0);
    // a[7] = std::complex<double>(1.0, 1.0);
    // a[8] = std::complex<double>(0.0, 0.0);
    //  Variables
    char jobvl = 'N';
    char jobvr = 'V';
    int n = dim;
    int lda = dim;
    int ldvl = dim;
    int ldvr = dim;
    int lwork_query = -1; // Query optimal workspace size
    int info;

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
        return 1;
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
        return 1;
    }

    // Print eigenvalues
    std::cout << "Eigenvalues:" << std::endl;
    for (int i = 0; i < dim; ++i)
    {
        std::cout << w[i] << std::endl;
    }

    return 0;
}
