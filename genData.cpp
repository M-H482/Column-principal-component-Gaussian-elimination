#include "helper.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <vector>
using namespace std;

void matVecMul(int n, float* A, float* x, float* b)
{
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        float sum = 0.0f;
        for (int j = 0; j < n; j++) {
            sum += REF(A, i, j, n) * x[j];
        }
        b[i] = sum;
    }
}

int main(int argc, char** argv)
{
    int N;
    char file[128];

    if (argc != 2) {
        printf("Usage: executable N\n");
        return 0;
    }

    N = atoi(argv[1]);

    sprintf(file, "./data/Axb_%d.txt", N);
    vector<float> A(N * N);
    vector<float> x(N);
    vector<float> b(N);

    int p = 10, q = 1;
    srand(time(NULL));
    for (int i = 0; i < N * N; i++)
        A[i] = 1.0 * (rand() % p + 1) / q;
    for (int i = 0; i < N; i++)
        x[i] = 1.0 * (rand() % p + 1) / q;

    matVecMul(N, A.data(), x.data(), b.data());

    FILE* fp = fopen(file, "w");
    fprintf(fp, "N = %d\n", N);

    for (int i = 0; i < N * N; i++)
        fprintf(fp, "%f\n", A[i]);
    for (int i = 0; i < N; i++)
        fprintf(fp, "%f\n", b[i]);
    fclose(fp);

    sprintf(file, "./data/Ref_%d.txt", N);
    fp = fopen(file, "w");

    for (int i = 0; i < N; i++)
        fprintf(fp, "%f\n", x[i]);
    fclose(fp);
    return 0;
}