// ukf.c
// Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

#include "ukf.h"

/*==============================================================================
                         STRUCTURE FUNCTION DEFINITIONS
==============================================================================*/
void ut(int dim, double* Z, double* W, double* mu, double** S, double a, double b, double k){

    // calculating weights
    int N = 2 * dim + 1; 
    double l = pow(a, 2.0) * (dim + k) - dim; 
    double c = 0.5 / (dim + l);
    W[0] = l / (dim + l);                          // mean weights
    W[N] = l / (dim + l) + (1 - pow(a, 2.0) + b);  // cov. weights
    for(int i = 1; i < N; i++){
        W[i] = c; W[i + N] = c; 
    }

    // calculating sigma points
    double** L = chol(dim, S);
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            L[i][j] = sqrt(dim + l) * L[i][j];
        }
        Z[i] = mu[i]; 
    }
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            Z[(i + 1) * dim + j]       = mu[j] + L[j][i]; 
            Z[(dim + i + 1) * dim + j] = mu[j] - L[j][i]; 
        }
    }

    // freeing dynamically allocated memory
    for(int i = 0; i < dim; i++){
        free(L[i]); 
    }
    free(L); 

    return;
}

/*==============================================================================
                         MISC. FUNCTION DEFINITIONS
==============================================================================*/
void update_moments_ukf(int dim, double* Z, double* W, double* mu, double** S){

    int N = 2 * dim + 1; 
    for(int j = 0; j < dim; j++){
        mu[j] = 0.0; 
        for(int i = 0; i < N; i++){
            mu[j] += W[i] * Z[i * dim + j];
        }
    }

    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            S[i][j] = 0.0;
        }
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < dim; j++){
            double diff_j = Z[i * dim + j] - mu[j];
            for(int k = 0; k < dim; k++){
                double diff_k = Z[i * dim + k] - mu[k];
                S[j][k] += W[i + N] * diff_j * diff_k;
            }
        }
    }

    return;
}

double norm(int dim, double* y){
    double sum = 0.0;
    for (int i = 0; i < dim; i++){
        sum += y[i] * y[i];
    }
    return sqrt(sum);
}

double calc_NED(int dim, double* Z, double* W, double* mu, double** S){

    // calculate mean and covariance of UKF
    update_moments_ukf(dim, Z, W, mu, S); 
    
    // calculate NED
    double* diff = init_v(dim); 
    for (int i = 0; i < dim; i++){
        diff[i] = Z[i] - mu[i];
    }
    double** L = chol(dim, S);
    if(L){
        double** Linv = inv(dim, L); 
        double* Linv_diff = mvtimes(dim, dim, Linv, diff); 
        double NED = norm(dim, Linv_diff);

        // cleanup
        for (int i = 0; i < dim; i++){
            free(L[i]);
            free(Linv[i]);
        }
        free(L); 
        free(Linv); 
        free(diff); 
        free(Linv_diff); 

        return NED;
    }else{
        return NAN; 
    }
}

void measurement_update_ukf(void (*h)(double t, double* x, double* y, double* coef), double** Z, double** W, double* y, double* coef, double** Sy, double t, int N, int dim, int y_dim, double a, double b, double k){
    
    double* mu = init_v(dim);
    double** S = init_m(dim, dim);

    update_moments_ukf(dim, *Z, *W, mu, S);

    // ------------------------------------------------------------
    // Transform sigma points through measurement model
    // ------------------------------------------------------------
    double** Y = init_m(N, y_dim);
    for(int i = 0; i < N; i++){
        h(t, &(*Z)[i * dim], Y[i], coef);
    }

    // ------------------------------------------------------------
    // Measurement mean
    // ------------------------------------------------------------
    double* y_mu = init_v(y_dim);
    for(int j = 0; j < y_dim; j++){

        if(j == 1){

            // Circular mean for theta
            double sin_sum = 0.0;
            double cos_sum = 0.0;

            for(int i = 0; i < N; i++){

                double theta = Y[i][j];

                sin_sum += (*W)[i] * sin(theta);
                cos_sum += (*W)[i] * cos(theta);
            }

            y_mu[j] = atan2(sin_sum, cos_sum);

        } else {

            // Ordinary mean
            y_mu[j] = 0.0;

            for(int i = 0; i < N; i++){
                y_mu[j] += (*W)[i] * Y[i][j];
            }
        }
    }

    // ------------------------------------------------------------
    // Measurement covariance and cross covariance
    // ------------------------------------------------------------
    double** P_yy = init_m(y_dim, y_dim);
    double** P_xy = init_m(dim, y_dim);
    for(int i = 0; i < y_dim; i++){
        for(int j = 0; j < y_dim; j++){
            P_yy[i][j] = 0.0;
        }
    }

    for(int i = 0; i < dim; i++){
        for(int j = 0; j < y_dim; j++){
            P_xy[i][j] = 0.0;
        }
    }

    for(int i = 0; i < N; i++){

        double dy[y_dim];

        // Wrapped measurement residuals
        for(int j = 0; j < y_dim; j++){
            dy[j] = measurement_residual(Y[i][j], y_mu[j], j);
        }

        // P_yy
        for(int j = 0; j < y_dim; j++){
            for(int l = 0; l < y_dim; l++){
                P_yy[j][l] += (*W)[i + N] * dy[j] * dy[l];
            }
        }

        // P_xy
        for(int j = 0; j < y_dim; j++){
            for(int l = 0; l < dim; l++){
                double dx = (*Z)[i * dim + l] - mu[l];
                P_xy[l][j] +=(*W)[i + N] * dx * dy[j];
            }
        }
    }

    // ------------------------------------------------------------
    // Add measurement noise
    // ------------------------------------------------------------
    for(int i = 0; i < y_dim; i++){
        for(int j = 0; j < y_dim; j++){
            P_yy[i][j] += Sy[i][j];
        }
    }

    // ------------------------------------------------------------
    // Kalman gain
    // ------------------------------------------------------------
    double** P_yy_inv = inv(y_dim, P_yy);
    double** K = mmtimes(dim, y_dim, y_dim, P_xy, P_yy_inv);

    // ------------------------------------------------------------
    // Measurement innovation
    // ------------------------------------------------------------
    double* diff_y = init_v(y_dim);
    for(int i = 0; i < y_dim; i++){
        diff_y[i] = measurement_residual(y[i], y_mu[i], i);
    }

    // ------------------------------------------------------------
    // State mean update
    // ------------------------------------------------------------
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < y_dim; j++){
            mu[i] += K[i][j] * diff_y[j];
        }
    }

    // ------------------------------------------------------------
    // Covariance update
    // ------------------------------------------------------------
    double** KT = transpose(dim, y_dim, K);
    double** K_Pyy = mmtimes(dim, y_dim, y_dim, K, P_yy);
    double** K_Pyy_KT = mmtimes(dim, y_dim, dim, K_Pyy, KT);
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            S[i][j] -= K_Pyy_KT[i][j];
        }
    }

    // ------------------------------------------------------------
    // Force symmetry
    // ------------------------------------------------------------
    for(int i = 0; i < dim; i++){
        for(int j = 0; j <= i; j++){
            double sym_val = 0.5 * (S[i][j] + S[j][i]);
            S[i][j] = sym_val;
            S[j][i] = sym_val;
        }
    }

    // ------------------------------------------------------------
    // Regenerate UKF sigma points
    // ------------------------------------------------------------
    ut(dim, *Z, *W, mu, S, a, b, k);

    // ------------------------------------------------------------
    // Cleanup
    // ------------------------------------------------------------
    free(mu);
    free(y_mu);
    free(diff_y);
    for(int i = 0; i < y_dim; i++){
        free(P_yy[i]);
        free(KT[i]);
        free(P_yy_inv[i]);
    }
    free(P_yy);
    free(KT);
    free(P_yy_inv);
    for(int i = 0; i < dim; i++){
        free(S[i]);
        free(P_xy[i]);
        free(K[i]);
        free(K_Pyy[i]);
        free(K_Pyy_KT[i]);
    }
    free(S);
    free(P_xy);
    free(K);
    free(K_Pyy);
    free(K_Pyy_KT);
    for(int i = 0; i < N; i++){
        free(Y[i]);
    }
    free(Y);

    return;
}
