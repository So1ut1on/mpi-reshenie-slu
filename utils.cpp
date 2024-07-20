#include "header.h"

void menu(char *argv[])
{
    printf("\nUsage: %s <args>\n\n", argv[0]);
    printf(" ___________________________Menu___________________________\n");
    printf("|  --help or -h - help menu                                |\n");
    printf("|  --print or -p - print matrix                            |\n");
    printf("|  -max_print=<int> - max amount of print (%d is default)  |\n",
           MAX_PRINT);
    printf("|  -eps=<double> - epsilon (default is 1e-16)              |\n");
    printf("|  -mode=<gen/read> - generate or read from file           |\n");
    printf("|  if mode is gen:                                         |\n");
    printf("|        -size=<int> - size of matr.Should be > 0          |\n");
    printf("|        -test_num=<0-2> - type of gen (%d is default)      |\n",
           TEST_NUM);
    printf("|  if mode is read:                                        |\n");
    printf("|        -path=[path/to/file] - path to file               |\n");
    printf(" __________________________________________________________\n");
}

double count_norm(double *a, int n, int my_rank, int p)
{
    double max = 0;
    double norm;
    for (int i = 0; i < n; i++) {
        int i_loc = i / p;
        if (i_loc == my_rank) {
            double s = 0;
            for (int j = 0; j < n; j++) {
                s += fabs(a[j + i_loc * n]);
            }
            if (s > max) {
                max = s;
            }
        }
    }
    MPI_Allreduce(&max, &norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return norm;
}

double count_error(double *r, double *row, int n, int my_rank, int p)
{
    double max = 0;
    double err;
    MPI_Status st;
    for (int i = 0; i < n; i++) {
        int i_loc = i / p;
        if (my_rank == 0) {
            if (my_rank == i % p) {
                row[i] = r[i_loc];
            } else {
                MPI_Recv(row + i, 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &st);
            }
        } else if (my_rank == i % p) {
            MPI_Send(r + i_loc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    for (int i = 0; i < n; i++) {
        int i_loc = i / p;
        if (i % p == my_rank) {
            double t = fabs(r[i_loc] - ((i + 1) % 2));
            if (t > max) {
                max = t;
            }
        }
    }
    MPI_Allreduce(&max, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return err;
}

double count_resid(const double *a, double *r, double *row, int n, int my_rank,
                   int p)
{
    double *rhs = new double[n];
    double res;
    double max = 0;
    MPI_Status st;
    for (int i = 0; i < n; i++) {
        int i_loc = i / p;
        if (my_rank == 0) {
            if (my_rank == i % p) {
                rhs[i] = row[i_loc];
            } else {
                MPI_Recv(rhs + i, 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &st);
            }
        } else if (my_rank == i % p) {
            MPI_Send(row + i_loc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Bcast(rhs, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i = 0; i < n; i++) {
        int i_loc = i / p;
        if (i % p == my_rank) {
            double s = 0;
            for (int j = 0; j < n; j++) {
                s += a[j + i_loc * n] * rhs[j];
            }

            double t = fabs(s - r[i_loc]);
            if (t > max) {
                max = t;
            }
        }
    }

    MPI_Allreduce(&max, &res, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    delete[] rhs;
    return res;
}

void shift_back(double *r, double *row, const int *transpos, int n, int my_rank,
                int p)
{
    MPI_Status st;
    for (int i = 0; i < n; i++) {
        int i_loc = i / p;
        if (my_rank == 0) {
            if (my_rank == i % p) {
                row[i] = r[i_loc];
            } else {
                MPI_Recv(row + i, 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &st);
            }
        } else if (my_rank == i % p) {
            MPI_Send(r + i_loc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    if (my_rank == 0) {
        for (int i = n - 1; i >= 0; i--) {
            if (i != transpos[i]) {
                double t = row[i];
                row[i] = row[transpos[i]];
                row[transpos[i]] = t;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        int i_loc = i / p;
        if (my_rank == 0) {
            if (my_rank == i % p) {
                r[i_loc] = row[i];
            } else {
                MPI_Send(row + i, 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD);
            }
        } else if (my_rank == i % p) {
            MPI_Recv(r + i_loc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &st);
        }
    }
}
