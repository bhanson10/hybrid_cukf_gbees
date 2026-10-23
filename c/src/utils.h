#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <omp.h>
#include <lapacke.h>
#include <cblas.h>
#include <gsl/gsl_cdf.h>

/*==============================================================================
                            FUNCTION DEFINITIONS
==============================================================================*/

double* init_v(int N);

int* init_v_i(int N);

double** init_m(int N, int M);

int** init_m_i(int N, int M);

double** cast_matrix(int N, int M, double matrix[N][M]);

double* mvtimes(int N, int M, double** matrix_a, const double* vector_b);

double** mmtimes(int N, int M, int L, double** matrix_a, double** matrix_b);

void print_vector(int N, double* vector);

void print_vector_T(int N, double* vector);

void print_matrix(int N, int M, double** matrix);

void print_matrix_i(int N, int M, int** matrix);

void read_measurement(int dim, const char* M_FILE, double* mu, double** S, double* T);

void read_measurement_p(int dim, const char* M_FILE, double* mu, double** S, double* JD_0, double* T);

void read_measurement_2_p(int dim, const char* M_FILE, double* mu, double* T);

void read_measurement_g(int dim, int y_dim, const char* M_FILE, double* x_g, double* mu, double** S, double* T);

void read_gbees(int N, int dim, const char* M_FILE, double** X, double* P);

char* JD_to_cal(double JD);

double dot(double* vector_a, double* vector_b, int size);

double** chol(int dim, double** S);

int chol2(int dim, double* S, double* L);

double** transpose(int N, int M, double** matrix);

void swap_rows(int dim, double** matrix, int row1, int row2);

void swap_vector_rows(int dim, double** vector, int row1, int row2);

double** inv(int dim, double** mat);

double generate_standard_normal(void);

double generate_standard_normal_p(unsigned int* seed);

double* mvnrnd(int N, int M, const double* mu, double** S);

double* mvnrnd_p(int N, int M, const double* mu, double** S, unsigned int* seed);

void RK4(void (*f)(double t, double* x, double* xp, double* coef), int dim, double* x, double* coef, double t, double dt);

double gauss_probability(int dim, double* x, double* mu, double** S);

double* mean(int N, int M, double* X, double* p);

double** cov(int N, int M, double* X, double* p);

double mahalanobis_distance(double *xi, double *xj, double **Sigma_inv, int dim);

double* cholesky_upper(int dim, double** S);

void solve_Rt(int dim, const double* R, const double* b, double* out);

double sq_euclid(int dim, const double* a, const double* b);

double* flatten_X(int N, int dim, double** X);

double calc_HZ_fast(int N, int dim, double* X, double* p, double* mu, double** S);

int logspace_fill(double a, double b, int n, double* out);

double* logspace(double a, double b, int n);

int count_txt_files(const char *path);

void measurement_update_gbees(void (*h)(double t, double* x, double* y, double* coef), double** X, double* P, double* y, double* coef, double** Sy, double t, double thresh, int N, int dim, int y_dim, double*** Xp, double** Pp, int* Np);

double wrap_angle(double a);

double measurement_residual(double a, double b, int index);

double norm_inf(double* g, int p);

double norm_frob(double* g, int p);

int expm(double* A, double* E, int n);

double** dt_process_noise(void (*J)(double t, double* x, double** F), double t, double** Q, double* x, int dim);

int parse_gbees(const char* filename, int dim, double** P, double*** X, int* N);

void weighted_mean_cov(double** X, double* P, int N, int dim, double* mu, double** Sigma);

void weighted_mean_cov_2(double* X, double* P, int N, int dim, double* mu, double** Sigma);

void write_pdf(char* FILE_NAME, double t, int N, int dim, double* X, double* P);
// void write_pdf(char* file_str, int N, int dim, double* X, double* P);

void write_pdf_2(char* file_str, int N, int dim, double* X, double* P, double* mu, double** S);

double det(int dim, const double* M);

void calc_b1_b2(int N, int dim, double* X, double* mu, double** S, double* b1, double* b2);

double calc_b1(int N, int dim, double* X, double* mu, double** S);

double calc_b2(int N, int dim, double* X, double* mu, double** S);

bool check_b1_b2(double b1, double b2, int N, int dim, double z1_s, double z2_s);

bool check_b1(double b1, int N, int dim, double z1_s);

bool check_b2(double b2, int N, int dim, double z2_s);

double frob_norm(int n, int m, double** A);

#endif // UTILS_H
