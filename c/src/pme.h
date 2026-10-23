#ifndef PME_H
#define PME_H

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include "cukf.h"
#include "utils.h"

/*==============================================================================
                              FUNCTION DEFINITIONS
==============================================================================*/
void dual_objective(const double* l, int p, int N, double* Phi, double* Mi, double w, double* f, double* g, double* H, double* s, double* e_scaled, double* Phi_scaled, bool grad, bool hess);

double comp_count(int r, int k);

void gen_E(int d, int M, double idx, int* vec);

void cut2mom(double** Z, double* W, int N, int M, int dim, int p, double* Mi, int** E);

double* solve_pme_eqns(double* Mi, int** E, double* X, double vol, int N, int dim, int p);

int whiten_cut(double** Z, double* W, int N, int M, int dim, double** Zp, double* mu, double* L);

double* phi(double* X, int N, double* lambda, int** E, int p, int dim);

void pme_cut2pdf(double** Z, double* W, int N, int CUT_M, int dim, int PME_M, int n_quad, int n_eval, double** X, double** P);

#endif // CUKF_H
