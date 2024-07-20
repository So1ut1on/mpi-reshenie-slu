#include "header.h"
double fun(int n, int i, int j, int test_num)
{
    ++i;
    ++j;
    switch (test_num) {
    case 1:
        return static_cast<double>(n - std::max(i, j) + 1);
    case 2:
        return static_cast<double>(std::max(i, j));
    case 3:
        return static_cast<double>(abs(i - j));
    case 4:
        return 1. / (i + j - 1);
#if 0
    case 0:
        return static_cast<double>(i == j);
    case 1:
        return fabs(i - j);
    case 2:
        return 1. / fabs(i + j + 1);
#endif
    }
    return i;
}

void generate_matrix(double *a, double *r, int n, int my_rank, int p,
                     int test_num)
{
    for (int i = my_rank; i < n; i += p) {
        int i_loc = i / p;
        r[i_loc] = 0;
        for (int j = 0; j < n; j++) {
            a[j + i_loc * n] = fun(n, i, j, test_num);
            if (j % 2 == 0) {
                r[i_loc] += a[j + i_loc * n];
            }
        }
    }
}

void generate_rhs_by_matrix(const double *a, double *r, int n, int my_rank,
                            int p)
{
    for (int i = my_rank; i < n; i += p) {
        int i_loc = i / p;
        r[i_loc] = 0;
        for (int j = 0; j < n; j++) {
            if (j % 2 == 0) {
                r[i_loc] += a[j + i_loc * n];
            }
        }
    }
}

int alloc_memory(double **a, double **r, double **row, int n, int my_rank,
                 int p)
{
    if (n > MAX_TO_ALLOC) {
        return -2;
    }
    int size = n / p;
    if (my_rank < n % p) {
        size++;
    }
    *a = new double[static_cast<size_t>(size) * static_cast<size_t>(n)];
    *r = new double[size];
    *row = new double[n + 1];

    if ((*a == nullptr) || (*r == nullptr) || (*row == nullptr)) {
        return -2;
    }
    return 0;
}

int read_matrix(double *a, double *row, int n, const char *file_name,
                int my_rank, int p)
{
    FILE *fp;
    int err = 0;
    MPI_Status st;
    if (my_rank == 0) {
        fp = fopen(file_name, "r");
        if (fp == nullptr) {
            err = -3;
        }
        MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (fp == nullptr) {
            return -3;
        }
        for (int i = 0; i < n; i++) {
            int i_loc = i / p;
            for (int j = 0; j < n; j++) {
                if (fscanf(fp, "%lf", row + j) != 1) {
                    // printf ("Wrong data or not enough elements\n");
                    err = -1;
                    MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    return -3;
                }
            }
            MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (i % p == 0) {
                memcpy(a + i_loc * n, row, n * sizeof(double));
                // r[i_loc] = row[n];
            } else {
                MPI_Send(row, n, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD);
                // MPI_Send(row + n, 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD);
            }
        }
        fclose(fp);
    } else {
        MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (err < 0) {
            return -3;
        }
        for (int i = 0; i < n; i++) {
            MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (err < 0) {
                return -3;
            }
            if (i % p == my_rank) {
                int i_loc = i / p;
                MPI_Recv(a + i_loc * n, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
                         &st);
                // MPI_Recv(r + i_loc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
                // &st);
            }
        }
    }

    return 0;
}

void print_matrix(double *a, double *row, int n, int my_rank, int p,
                  int max_print, double *r)
{
    int max_n = (n < max_print) ? n : max_print;
    MPI_Status st;
    for (int i = 0; i < max_n; i++) {
        int i_loc = i / p;
        if (my_rank == 0) {
            if (i % p == 0) {
                for (int j = 0; j < max_n; j++) {
                    printf(" %.3le", a[j + i_loc * n]);
                }
                if (r != nullptr) {
                    printf(" %.3le", r[i_loc]);
                }
                printf("\n");
            } else {
                MPI_Recv(row, max_n, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &st);
                for (int j = 0; j < max_n; j++) {
                    printf(" %.3le", row[j]);
                }
                if (r != nullptr) {
                    double var;
                    MPI_Recv(&var, 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD,
                             &st);
                    printf(" %.3le", var);
                }
                printf("\n");
            }
        } else if (my_rank == i % p) {
            MPI_Send(a + i_loc * n, max_n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            if (r != nullptr) {
                MPI_Send(r + i_loc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }
    }
}

void print_vec(double *r, int n, int my_rank, int p, int max_print)
{
#if 0
  if (my_rank == 0)
    printf ("Vector\n");
#endif
    int max_n = (n < max_print) ? n : max_print;
    MPI_Status st;
    for (int i = 0; i < max_n; i++) {
        int i_loc = i / p;
        if (my_rank == 0) {
            if (i % p == 0) {
                printf("%e ", r[i_loc]);
            } else {
                double var;
                MPI_Recv(&var, 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &st);
                printf("%e ", var);
            }
        } else if (my_rank == i % p) {
            MPI_Send(r + i_loc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    if (my_rank == 0) {
        printf("\n");
    }
}
