// pf.c
// Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

#include "utils.h"

double det_2x2(double** S){
    return S[0][0] * S[1][1] - S[0][1] * S[1][0];
}

double det_3x3(double** S){
    return S[0][0] * (S[1][1] * S[2][2] - S[1][2] * S[2][1])
         - S[0][1] * (S[1][0] * S[2][2] - S[1][2] * S[2][0])
         + S[0][2] * (S[1][0] * S[2][1] - S[1][1] * S[2][0]);
}

double log_gauss_probability(int dim, double* x, double* mu, double** S, double** S_inv){ 
    double diff[dim];
    for(int i=0; i<dim; i++) diff[i] = x[i] - mu[i];

    double* S_inv_x = mvtimes(dim, dim, S_inv, diff);
    double dot_prod = dot(diff, S_inv_x, dim);
    double det_S; 
    if(dim == 2){
        det_S = det_2x2(S); 
    }else if(dim == 3){
        det_S = det_3x3(S); 
    }else{
        fprintf(stderr, "log_gauss_probability: unavailable measurement dimension for determinant.");
        exit(EXIT_FAILURE);
    }
    free(S_inv_x);

    return -0.5 * (dim * log(2 * M_PI) + log(det_S) + dot_prod);
}

int find_index_binary(double* cdf, int N, double u){

    int lo = 0, hi = N-1;
    while(lo < hi){
        int mid = lo + (hi - lo)/2;
        if(cdf[mid] >= u) hi = mid;
        else lo = mid + 1;
    }
    return lo;
}

void measurement_update_pf(void (*h)(double t, double* x, double* y, double* coef), double* X, double* P, double* y, double* coef, double** Sy, double t, double ess_t, int N, int dim, int y_dim){

    double* l = init_v(N);
    double* log_l = init_v(N);
    double* y_i = init_v(y_dim); 

    // 1. Compute LOG-likelihoods
    double** Sy_inv = inv(y_dim, Sy);
    double max_log_l = -DBL_MAX;
    for(int i = 0; i < N; i++){
        h(t, &X[i * dim], y_i, coef); 
        // Return log-probability instead of raw probability
        log_l[i] = log_gauss_probability(y_dim, y_i, y, Sy, Sy_inv); 
        if(log_l[i] > max_log_l){
            max_log_l = log_l[i]; // Track maximum log-likelihood for stability
        }
    }

    // 2. Compute normalized weights P[i] using Log-Sum-Exp
    double sum_exp = 0.0;
    for(int i = 0; i < N; i++){
        // Subtract max_log_l to prevent underflow/overflow
        l[i] = exp(log_l[i] - max_log_l); 
        sum_exp += l[i];
    }

    // 3. Compute Normalized Weights and ESS
    double sum_sq_P = 0.0; 
    for(int i = 0; i < N; i++){
        P[i] = l[i] / sum_exp; 
        sum_sq_P += P[i] * P[i]; 
    }

    double ess = 1.0 / sum_sq_P; 
    printf("ESS: %.16f\n", ess);

    // // effective sample size
    // double ess = 0.0; 
    // for(int i = 0; i < N; i++){
    // P[i] = l[i] / sum; 
    // ess += pow(P[i], 2.0); 
    // }
    // ess = 1.0 / ess; 
    // printf("ESS: %.16f\n", ess); 

    if(ess < ess_t){

        // cdf
        double* cdf = init_v(N);
        cdf[0] = P[0];
        for(int i = 1; i < N; i++){
            cdf[i] = cdf[i - 1] + P[i]; 
        }
        cdf[N - 1] = 1.0; 

        int* re_idx = init_v_i(N);
        double** Xp = init_m(N, dim);
        for(int i = 0; i < N; i++){
            double u = (double)rand() / RAND_MAX;
            int idx = find_index_binary(cdf, N, u);
            re_idx[i] = idx;
        }

        for(int i = 0; i < N; i++){
            for(int j = 0; j < dim; j++){
            Xp[i][j] = X[re_idx[i] * dim + j]; 
            }
        }

        for(int i = 0; i < N; i++){
            for(int j = 0; j < dim; j++){
                X[i * dim + j] = Xp[i][j]; 
            }
            P[i] = 1.0 / N; 
        }

        free(cdf);
        free(re_idx); 
        for(int i = 0; i < N; i++){
            free(Xp[i]); 
        }
        free(Xp); 
    }

    for(int i = 0; i < y_dim; i++){
        free(Sy_inv[i]); 
    }
    free(Sy_inv); 
    free(l); 
    free(log_l); 
    free(y_i); 

    return;
}
