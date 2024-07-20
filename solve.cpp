#include "header.h"

int solve(double *a, double *r, double *row, int *transpos, int n, int my_rank,
          int p, double eps, double norm)
{
    int size = n / p;
    if (my_rank < n % p) {
        size++;
    }
    // 1 step of gauss
    for (int step = 0; step < n; step++) {
        int sw;
        int i_loc = step / p;
        if (my_rank == step % p) {
            double max = fabs(a[step + i_loc * n]);
            sw = step;
            for (int j = step + 1; j < n; j++) {
                double t = fabs(a[j + i_loc * n]);
                if (t > max) {
                    max = t;
                    sw = j;
                }
            }
            if (max < eps * norm) {
                sw = -1;
            } else {
                double coeff = a[sw + i_loc * n];

                for (int j = step; j < n; j++) {
                    a[j + i_loc * n] /= coeff;
                }
                r[i_loc] /= coeff;
            }
        }
        MPI_Bcast(&sw, 1, MPI_INT, step % p, MPI_COMM_WORLD);
        if (sw < 0) {
            return -1;
        }
        transpos[step] = sw;

        if (sw != step) {
            for (int i = 0; i < size; i++) {
                double t = a[sw + i * n];
                a[sw + i * n] = a[step + i * n];
                a[step + i * n] = t;
            }
        }

        if (my_rank == step % p) {
            memcpy(row + step + 1, a + i_loc * n + step + 1,
                   (n - step - 1) * sizeof(double));
            row[n] = r[i_loc];
        }

        MPI_Bcast(row + step + 1, n - step, MPI_DOUBLE, step % p,
                  MPI_COMM_WORLD);
        int add = 0;
        if (step % p >= my_rank) {
            add = 1;
        }
        for (int i = step / p + add; i < size; i++) {
            double coeff = a[step + i * n];
            for (int j = step + 1; j < n; j++) {
                a[j + i * n] -= coeff * row[j];
            }

            r[i] -= coeff * row[n];
        }
    }

    // 2 step of gauss
    for (int step = n - 1; step > 0; step--) {
        double coeff = 0;
        if (step % p == my_rank) {
            coeff = r[step / p];
        }
        MPI_Bcast(&coeff, 1, MPI_DOUBLE, step % p, MPI_COMM_WORLD);
        int add = 0;
        if ((step - 1) % p < my_rank) {
            add = 1;
        }
        for (int i = (step - 1) / p - add; i >= 0; i--) {
            r[i] -= coeff * a[step + i * n];
        }
    }

    return 0;
}
