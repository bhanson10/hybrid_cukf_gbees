// ukf.h
// Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

#ifndef UKF_H
#define UKF_H

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
void ut(int dim, double* Z, double* W, double* mu, double** S, double a, double b, double k);

void update_moments_ukf(int dim, double* Z, double* W, double* mu, double** S);

double norm(int dim, double* y);

double calc_NED(int dim, double* Z, double* W, double* mu, double** S);

void measurement_update_ukf(void (*h)(double t, double* x, double* y, double* coef), double** Z, double** W, double* y, double* coef, double** Sy, double t, int N, int dim, int y_dim, double a, double b, double k);

#endif // UKF_H
