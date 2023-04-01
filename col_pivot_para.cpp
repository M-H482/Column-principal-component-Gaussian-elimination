#include "helper.h"
#include "mpi.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <sys/time.h>
using namespace std;

void printX(int rank, int nprcs, int n, float* x_sub);
void printY(int rank, int nprcs, int n, float* y_sub);
void printLU(int rank, int nprcs, int n, float* A_sub, float* b);

int find_max(int n, float* A, int j, int icol)
{
    int r = j;
    float v = fabsf(REF(A, j, icol, n));

    for (int i = j + 1; i < n; i++) {
        float u = fabsf(REF(A, i, icol, n));
        if (u > v) {
            v = u;
            r = i;
        }
    }
    return r;
}

void myswap(int n, int m, float* A, int i, int j)
{
    for (int k = 0; k < m; k++) {
        swap(REF(A, i, k, n), REF(A, j, k, n));
    }
}

void guass_elimintation(int rank, int nprcs, int n, float* A_sub, float* b)
{
    int icol = 0, m = n / nprcs;
    int rankp1 = (rank + 1) % nprcs, rankm1 = (rank - 1 + nprcs) % nprcs;
    float* f = new float[n];
    int* l = new int[n];

    for (int j = 0; j < n - 1; j++) {
        if (rank == j % nprcs) {
            l[j] = find_max(n, A_sub, j, icol);
            if (l[j] != j) {
                myswap(n, m, A_sub, l[j], j);
            }
            float pivot = REF(A_sub, j, icol, n);

            if (fabsf(pivot) < 1e-6) {
                printf("A is singular! exit now!\n");
                MPI_Abort(MPI_COMM_WORLD, -1);
            }

            for (int i = j + 1; i < n; i++) {
                REF(A_sub, i, icol, n) = REF(A_sub, i, icol, n) / pivot;
                f[i] = REF(A_sub, i, icol, n);
            }

            icol += 1;

            MPI_Send(&l[j], 1, MPI_INT, rankp1, 2323, MPI_COMM_WORLD);
            MPI_Send(f + j + 1, n - j - 1, MPI_FLOAT, rankp1, 8848, MPI_COMM_WORLD);
        } else {
            MPI_Status sta;

            MPI_Recv(&l[j], 1, MPI_INT, rankm1, 2323, MPI_COMM_WORLD, &sta);
            MPI_Recv(f + j + 1, n - j - 1, MPI_FLOAT,
                rankm1, 8848, MPI_COMM_WORLD, &sta);

            if (rankp1 != j % nprcs) {
                MPI_Send(&l[j], 1, MPI_INT, rankp1, 2323, MPI_COMM_WORLD);
                MPI_Send(f + j + 1, n - j - 1, MPI_FLOAT, rankp1, 8848, MPI_COMM_WORLD);
            }

            if (l[j] != j) {
                myswap(n, m, A_sub, l[j], j);
            }
        }

        for (int k = icol; k < m; k++) {
            for (int i = j + 1; i < n; i++) {
                REF(A_sub, i, k, n) -= f[i] * REF(A_sub, j, k, n);
            }
        }
    }

    if (rank == 0) {
        for (int i = 0; i < n - 1; i++) {
            swap(b[i], b[l[i]]);
        }
    }
    delete[] f;
    delete[] l;
}

void getY(int rank, int p, int n, float* L, float* b, float* y)
{ // Ly=b
    MPI_Status sta;
    float* u = new float[n];
    float* v = new float[n];
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            u[i] = b[i];
            v[i] = 0.0f;
        }
    } else {
        for (int i = 0; i < n; i++) {
            u[i] = 0.0f;
            v[i] = 0.0f;
        }
    }

    int k = 0, m = n / p;

    for (int i = rank; i < n; i += p) {
        if (i > 0) {
            MPI_Recv(v, p - 1, MPI_FLOAT,
                (i - 1) % p, 777, MPI_COMM_WORLD, &sta);
        }

        // y[k] = (u[i] + v[0]) / REF(L, i, k, m);
        y[k] = u[i] + v[0];

        for (int j = 0; j < p - 2 && i + 1 + j < n; j++) {
            v[j] = v[j + 1] + u[i + 1 + j] - REF(L, i + 1 + j, k, n) * y[k];
        }
        if (i + p - 1 < n) {
            v[p - 2] = u[i + p - 1] - REF(L, i + p - 1, k, n) * y[k];
        }

        if (i < n)
            MPI_Send(v, p - 1, MPI_FLOAT, (i + 1) % p, 777, MPI_COMM_WORLD);

        for (int j = i + p; j < n; j++) {
            u[j] = u[j] - REF(L, j, k, n) * y[k];
        }

        k = k + 1;
    }
    delete[] u;
    delete[] v;
}

void getX(int rank, int p, int n, float* U, float* y, float* x)
{ // Ux=y
    MPI_Status sta;
    float* u = new float[n];
    float* v = new float[n];

    for (int i = 0; i < n; i++) {
        v[i] = 0.0f;
    }

    for (int k = 0; k < n / p; k++) {
        MPI_Gather(y + k, 1, MPI_FLOAT, u + k * p, 1, MPI_FLOAT, p - 1, MPI_COMM_WORLD);
    }

    if (rank != p - 1) {
        for (int i = 0; i < n; i++) {
            u[i] = 0.0f;
        }
    }

    int m = n / p;
    int k = m - 1;

    for (int i = n - (p - rank); i >= 0; i -= p) {

        if (i < n - 1) {
            MPI_Recv(v, p - 1, MPI_FLOAT,
                (i + 1) % p, 888, MPI_COMM_WORLD, &sta);
        }

        x[k] = (u[i] + v[0]) / REF(U, i, k, n);

        for (int j = 0; j < p - 2 && i - 1 - j >= 0; j++) {
            v[j] = v[j + 1] + u[i - 1 - j] - REF(U, i - 1 - j, k, n) * x[k];
        }
        if (i - p + 1 >= 0) {
            v[p - 2] = u[i - p + 1] - REF(U, i - p + 1, k, n) * x[k];
        }

        if (i > 0) {
            MPI_Send(v, p - 1, MPI_FLOAT, (i - 1) % p, 888, MPI_COMM_WORLD);
        }

        for (int j = i - p; j >= 0; j--) {
            u[j] = u[j] - REF(U, j, k, n) * x[k];
        }

        k = k - 1;
    }
    delete[] u;
    delete[] v;
}

void scatterA(int rank, int nprcs, int n, float* A, float* A_sub)
{
#ifdef COL_MAJOR
    for (int k = 0; k < n / nprcs; k++) {
        MPI_Scatter(A + k * n * nprcs, n, MPI_FLOAT,
            A_sub + k * n, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
#endif
}

int main(int argc, char** argv)
{
    int rank, nprcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprcs);

    char file[128];
    int N;
    float *A, *x, *x_ref, *b;
    float *A_sub, *x_sub, *y_sub;

    A = x = x_ref = b = nullptr;
    A_sub = x_sub = y_sub = nullptr;

    if (argc != 2) {
        if (rank == 0)
            printf("Usage: executable file\n");
        MPI_Finalize();
        return 0;
    }

    if (rank == 0) {
        strcpy(file, argv[1]);
        init_element(file, N, A, x, b, x_ref);
        // for (int i = 0; i < N; i++) {
        //     for (int j = 0; j < N; j++) {
        //         printf("%f ", REF(A, i, j, N));
        //     }
        //     printf("\n");
        // }
        // printf("\n");
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    A_sub = new float[N * N / nprcs];
    x_sub = new float[N / nprcs];
    y_sub = new float[N / nprcs];

    scatterA(rank, nprcs, N, A, A_sub);

    double start, end, elapsed_time;
    start = MPI_Wtime();

    guass_elimintation(rank, nprcs, N, A_sub, b);

    getY(rank, nprcs, N, A_sub, b, y_sub);

    getX(rank, nprcs, N, A_sub, y_sub, x_sub);

    for (int k = 0; k < N / nprcs; k++) {
        MPI_Gather(x_sub + k, 1, MPI_FLOAT, x + k * nprcs, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    end = MPI_Wtime();
    elapsed_time = end - start;

    if (rank == 0) {
        printf("Elapsed time is %f us, %lf GFLOPS\n", elapsed_time * 1e6,
            (2.0 * N * N * N / 3.0 + 4.0 * N / 3.0) / elapsed_time / 1e9);
        if (check(N, x, x_ref, 1e-1)) {
            printf("The answer is right!\n");
        } else {
            printf("The answer is wrong!\n");
        }
    }

    delete[] A;
    delete[] x;
    delete[] x_ref;
    delete[] b;
    delete[] x_sub;
    delete[] y_sub;

    MPI_Finalize();
    return 0;
}

void printLU(int rank, int nprcs, int n, float* A_sub, float* b)
{
    float* A = nullptr;
    if (rank == 0) {
        A = new float[n * n];
    }

    for (int k = 0; k < n / nprcs; k++) {
        MPI_Gather(A_sub + k * n, n, MPI_FLOAT,
            &REF(A, 0, k * nprcs, n), n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        printf("after guass, [A b]\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%9f ", REF(A, i, j, n));
            }
            printf("%9f\n", b[i]);
        }
        printf("\n");
        delete[] A;
    }
}

void printY(int rank, int nprcs, int n, float* y_sub)
{
    float* y = nullptr;
    if (rank == 0) {
        y = new float[n];
    }

    for (int k = 0; k < n / nprcs; k++) {
        MPI_Gather(y_sub + k, 1, MPI_FLOAT, y + k * nprcs, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        printf("Y\n");
        for (int i = 0; i < n; i++) {
            printf("%f ", y[i]);
        }
        printf("\n\n");
        delete[] y;
    }
}

void printX(int rank, int nprcs, int n, float* x_sub)
{
    float* x = nullptr;
    if (rank == 0) {
        x = new float[n];
    }

    for (int k = 0; k < n / nprcs; k++) {
        MPI_Gather(x_sub + k, 1, MPI_FLOAT, x + k * nprcs, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        printf("X\n");
        for (int i = 0; i < n; i++) {
            printf("%f ", x[i]);
        }
        printf("\n\n");
        delete[] x;
    }
}
