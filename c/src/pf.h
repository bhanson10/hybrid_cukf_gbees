// pf.h
// Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

#ifndef PF_H
#define PF_H

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
double det_2x2(double** S);

double det_3x3(double** S);

double log_gauss_probability(int dim, double* x, double* mu, double** S, double** S_inv);

int find_index_binary(double* cdf, int N, double u);

void measurement_update_pf(void (*h)(double t, double* x, double* y, double* coef), double* X, double* P, double* y, double* coef, double** Sy, double t, double ess_t, int N, int dim, int y_dim);

#endif // PF_H
