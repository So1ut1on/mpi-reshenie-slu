#ifndef HEADER_H
#define HEADER_H
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include "mpi.h"
#include <memory>
#include <math.h>
#include <string>
#include <sys/time.h>
#include <sys/resource.h>
using namespace std;

double fun (int n, int i, int j, int test_num);
void generate_matrix (double *a, double *r, int n, int my_rank, int p, int test_num);
void generate_rhs_by_matrix (const double *a, double *r, int n, int my_rank, int p);
int alloc_memory (double **a, double **r, double **row, int n, int my_rank, int p);
int read_matrix (double *a, double *row, int n, const char *file_name, int my_rank, int p);
void print_matrix (double *a, double *row, int n, int my_rank, int p, int max_print, double *r = nullptr);
void print_vec (double *r, int n, int my_rank, int p, int max_print);
//--------------------------------------------------
//--------------------------------------------------
int solve (double *a, double *r, double *row, int *transpos, int n, int my_rank, int p, double eps, double norm);
double count_norm (double *a, int n, int my_rank, int p);
double count_resid (const double *a, double *r, double *row, int n, int my_rank, int p);
double count_error (double *r, double *row, int n, int my_rank, int p);
void shift_back (double *r, double *row, const int *transpos, int n, int my_rank, int p);
//--------------------------------------------------
//--------------------------------------------------
void menu (char *argv[]);
//--------------------------------------------------
//--------------------------------------------------
#define MAX_PRINT 10
#define TEST_NUM 1
#define MAX_TEST 4
#define MAX_TO_ALLOC 20000
#endif // HEADER_H
