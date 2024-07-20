#include "header.h"

double get_earth_time()
{
    const double precision = 1e-6;
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * precision;
}

double get_process_time()
{
    const double precision = 1e-6;
    struct rusage time;
    getrusage(RUSAGE_SELF, &time);
    return (double)time.ru_utime.tv_sec +
           (double)time.ru_utime.tv_usec * precision;
}

int main(int argc, char *argv[])
{

    int my_rank;
    int p;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    double *a = nullptr;
    double *r = nullptr;
    double *row = nullptr;
    const char *name = nullptr;
    int n = 0;
    int test = TEST_NUM;
    const double eps = 1e-16;
    double norm;
    int max_pr;
    bool print = false;

    const int min_args = 4;
    const int max_args = 5;
    if (argc < min_args || argc > max_args) {
        MPI_Finalize();
        return -1;
    }

    if ((n = atoi(argv[1])) <= 0 || (max_pr = atoi(argv[2])) <= 0 ||
        (test = atoi(argv[3])) < 0 || test > MAX_TEST ||
        (argc == 4 && test == 0) || max_pr > n) {
        MPI_Finalize();
        return -1;
    }

    if (argc == max_args && test != 0) {
        MPI_Finalize();
        return -1;
    }

    string mode;
    if (argc == max_args) {
        name = argv[4];
        mode = "read";
    } else {
        mode = "gen";
    }

    print = max_pr != 0;

    int err = 0;
    err = alloc_memory(&a, &r, &row, n, my_rank, p);
    if (err < 0) {
        delete[] a;
        delete[] r;
        delete[] row;

        MPI_Finalize();
        return err;
    }
    if (my_rank == 0) {
        // printf("allocation done\n");
    }
    if (mode == "gen") {
        generate_matrix(a, r, n, my_rank, p, test);
    }
    if (mode == "read") {
        err = read_matrix(a, row, n, name, my_rank, p);
        generate_rhs_by_matrix(a, r, n, my_rank, p);
    }
    if (my_rank == 0) {
        // printf("reading done: %d\n", err);
    }
    if (err == 0) {
        if (print) {
            print_matrix(a, row, n, my_rank, p, max_pr, r);
            // print_vec (r, n, my_rank, p, max_pr);
        }
        norm = count_norm(a, n, my_rank, p);
        int *transpos = new int[n];

        double global_time = get_earth_time();
        double time = get_process_time(); // MPI_Wtime ();
        int res = solve(a, r, row, transpos, n, my_rank, p, eps, norm);
        time = get_process_time() - time; // MPI_Wtime () - time;
        global_time = get_earth_time() - global_time;

        if (res < 0) {
            err = -4;
            // printf("solver: %d, err = %d\n", res, err);
#if 0
          if (my_rank == 0)
            {
              printf ("Det = 0\n");
            }
#endif
        } else {
            shift_back(r, row, transpos, n, my_rank, p);
            if (my_rank == 0) {
                printf("Solution: ");
            }
            print_vec(r, n, my_rank, p, max_pr);

            int size_n = n / p;
            if (my_rank < n % p) {
                size_n++;
            }
            double *reserve = new double[size_n];
            memcpy(reserve, r, size_n * sizeof(double));

            // Print process time
            MPI_Status st;
            for (int proc_num = 0; proc_num < p; ++proc_num) {
                double proc_time;
                if (my_rank == 0) {
                    if (proc_num == 0) {
                        printf("Time: %lf\n", global_time);
                        proc_time = time;
                    } else {
                        MPI_Recv(&proc_time, 1, MPI_DOUBLE, proc_num, 0,
                                 MPI_COMM_WORLD, &st);
                    }
                    printf("Time of process %d: %lf\n", proc_num, proc_time);
                } else if (my_rank == proc_num) {
                    MPI_Send(&time, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                }
            }

            // Print Error
            // if (mode == "gen")
            {
                double error = count_error(r, row, n, my_rank, p);
                generate_matrix(a, r, n, my_rank, p, test);
                if (my_rank == 0) {
                    printf("Error: %e\n", error);
                }
            }
            if (mode == "read") {
                read_matrix(a, row, n, name, my_rank, p);
                generate_rhs_by_matrix(a, r, n, my_rank, p);
            }
            double res = count_resid(a, r, reserve, n, my_rank, p);

            // printf ("Time of process %d: %lf\n", my_rank, time);
            if (my_rank == 0) {
                printf("Residual: %e\n", res);
            }
            delete[] reserve;
        }
        delete[] transpos;
    } else {
#if 0
      if (my_rank == 0)
        {
          printf ("Problems with reading\n");
        }
#endif
    }
    delete[] a;
    delete[] r;
    delete[] row;
    MPI_Finalize();
    return err;
}
