#ifndef UTILS_ASSIST_H
#define UTILS_ASSIST_H

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
#include "rebound.h"
#include "utils.h"
#include "ukf.h"
#include "cukf.h"

/*==============================================================================
                            FUNCTION DEFINITIONS
==============================================================================*/

double* mean_particles(int N, int M, const struct reb_particle* particles, const double* p);

double** cov_particles(int N, int M, const struct reb_particle* particles, const double* p, const double* mu);

double calc_HZ_fast_particles(int N, int dim, const struct reb_particle* particles, double* p, double* mu, double** S);

double calc_b1_particles(int N, int dim, const struct reb_particle* particles, double* mu, double** S);

double calc_b2_particles(int N, int dim, const struct reb_particle* particles, double* mu, double** S);

void update_moments_ukf_particles(int dim, const struct reb_particle* particles, double* W, double* mu, double** S);

double calc_NED_particles(int dim, const struct reb_particle* particles, double* W, double* mu, double** S);

double* reconstruct_cut_particles(int dim, int M1, int M2, const struct reb_particle* particles, const double* W);

double calc_bt1_particles(int N, int dim, int M, const struct reb_particle* particles, double* W, double* mu, double** S);

double calc_bt2_particles(int N, int dim, int M, const struct reb_particle* particles, double* W, double* mu, double** S);

#endif // UTILS_ASSIST_H