// ukf.h
// Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

#ifndef CUKF_H
#define CUKF_H

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include "utils.h"

/*==============================================================================
                              FUNCTION DEFINITIONS
==============================================================================*/
double** eye(int dim);

double** sqrtm(int dim, double** P);

void cut4(int dim, double* Z, double* W, double* mu, double** S);

void cut6(int dim, double* Z, double* W, double* mu, double** S);

double** principal_points(int dim);

double** conj_N_points(int dim);

double** conj_M_points(int dim, int order);

int nchoosek(int N, int M);

double* reconstruct_cut(int dim, int M1, int M2, double* Z, double* W);

double* outer_M(double* v, int d, int M);

void measurement_update_cukf(void (*h)(double t, double* x, double* y, double* coef), double** Z, double** W, double* y, double* coef, double** Sy, double t, int N, int dim, int y_dim, int M);

void measurement_update_cukf_pme(void (*h)(double t, double* x, double* y, double* coef), double** Z, double** W, double* X, double* P, int num_pts, double* y, double* coef, double** Sy, double t, int N, int dim, int y_dim, int M);

void calc_bt1_bt2(int N, int dim, int M, double* Z, double* W, double* mu, double** S, double* bt1, double* bt2);

double calc_bt1(int N, int dim, int M, double* Z, double* W, double* mu, double** S);

double calc_bt2(int N, int dim, int M, double* Z, double* W, double* mu, double** S);

#endif // CUKF_H
