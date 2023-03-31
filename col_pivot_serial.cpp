#include "helper.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <sys/time.h>
using namespace std;

int find_max(int n, float* A, int j)
{
    int r = j;
    float v = fabsf(REF(A, j, j, n));
    for (int i = j; i < n; i++) {
        float u = fabsf(REF(A, i, j, n));
        if (u > v) {
            v = u;
            r = i;
        }
    }
    return r;
}

void myswap(int n, float* A, float* b, int i, int j)
{
    for (int k = 0; k < n; k++) {
        swap(REF(A, i, k, n), REF(A, j, k, n));
    }
    swap(b[i], b[j]);
}

void guass_elimination(int n, float* A, float* b)
{
    for (int j = 0; j < n - 1; j++) {
        int l = find_max(n, A, j);
        if (l != j)
            myswap(n, A, b, l, j);
        if (REF(A, j, j, n) > -1e-6 && REF(A, j, j, n) < 1e-6) {
            printf("A is singular! exit now!\n");
            return;
        }

        float Ajj = REF(A, j, j, n);

        for (int i = j + 1; i < n; i++) {
            REF(A, i, j, n) /= Ajj;
        }

        for (int k = j + 1; k < n; k++) {
            float Ajk = REF(A, j, k, n);

            for (int i = j + 1; i < n; i++) {
                REF(A, i, k, n) -= REF(A, i, j, n) * Ajk;
            }
        }
    }
}

void getY(int n, float* A, float* b, float* y)
{ // LY=b
    for (int j = 0; j < n; j++) {
        float yj = y[j] = b[j];

        for (int i = j + 1; i < n; i++) {
            b[i] = b[i] - REF(A, i, j, n) * yj;
        }
    }
}

void getX(int n, float* A, float* y, float* x)
{ // Ux=y
    for (int j = n - 1; j >= 0; j--) {
        x[j] = y[j] / REF(A, j, j, n);
        float xj = x[j];

        for (int i = 0; i < j; i++) {
            y[i] = y[i] - REF(A, i, j, n) * xj;
        }
    }
}

int main(int argc, char** argv)
{

    char file[128] = "./data/Axb_4.txt";
    int N;
    float *A, *x, *x_ref, *b, *y;

    A = x = x_ref = b = y = nullptr;

    if (argc != 2) {
        printf("Usage: executable file\n");
        return 0;
    }
    strcpy(file, argv[1]);

    init_element(file, N, A, x, b, x_ref);
    y = new float[N];

    double elapsed_time;
    timeval start, end;
    gettimeofday(&start, nullptr);

    guass_elimination(N, A, b);

    getY(N, A, b, y);

    getX(N, A, y, x);

    gettimeofday(&end, nullptr);
    elapsed_time = (end.tv_sec - start.tv_sec) * 1e9 + (end.tv_usec - start.tv_usec);

    printf("Elapsed time is %lf us, %lf GFLOPS\n", elapsed_time * 1e6,
        (2.0 * N * N * N / 3.0 + 4.0 * N / 3.0) / elapsed_time / 1e9);

    if (check(N, x, x_ref, 1e-1)) {
        printf("The answer is right!\n");
    } else {
        printf("The answer is wrong!\n");
    }
    // sprintf(file, "./data/Ans_%d.txt", N);
    // FILE* fp = fopen(file, "w");

    // for (int i = 0; i < N; i++) {
    //     fprintf(fp, "%f\n", x[i]);
    // }

    // fclose(fp);

    delete[] A;
    delete[] x;
    delete[] x_ref;
    delete[] b;
    delete[] y;
    return 0;
}