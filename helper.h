#ifndef __HELPER_H__
#define __HELPER_H__

#include <cmath>
#include <iostream>

#define COL_MAJOR

#ifdef COL_MAJOR
#define REF(A, i, j, ld) A[(j) * (ld) + (i)]
#endif

#ifdef ROW_MAJOR
#define REF(A, i, j, ld) A[(i) * (ld) + (j)]
#endif

void init_element(char* file, int& N, float*& A, float*& x, float*& b, float*& x_ref)
{
    FILE* fp = fopen(file, "r");
    int a, r;

    r = fscanf(fp, "N = %d\n", &a);
    N = a;

    printf("N = %d\n", N);
    A = new float[N * N], x = new float[N], b = new float[N];
    x_ref = new float[N];

    for (int i = 0; i < N * N; i++)
        r = fscanf(fp, "%f", &A[i]);

    for (int i = 0; i < N; i++) {
        r = fscanf(fp, "%f", &b[i]);
    }

    fclose(fp);
}

bool check(int n, float* ans, float* ref, float threshold)
{
    char file[128];
    sprintf(file, "./data/Ref_%d.txt", n);

    FILE* fp = fopen(file, "r");

    for (int i = 0; i < n; i++) {
        int ret = fscanf(fp, "%f", &ref[i]);
    }
    fclose(fp);

    int j = 0;
    for (int i = 0; i < n; i++) {
        if (isnanf(ans[i])) {
            printf("check error: ans[%d] is nan!\n", i);
            return 1;
        }
        if (fabsf(ans[i] - ref[i]) > threshold) {
            printf("Error on index %d, ans = %f and ref = %f\n", i, ans[i], ref[i]);
            if (++j == 3) {
                printf("...\n");
                break;
            }
        }
    }
    return j == 0;
}

#endif