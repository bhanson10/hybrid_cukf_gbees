// utils_assist.c
// Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

#include "utils_assist.h"

double* mean_particles(int N, int M, const struct reb_particle* particles, const double* p){
    if (N <= 0 || M <= 0 || !particles) return NULL;

    double* m = init_v(M); 
    if (!m) return NULL;

    const double default_weight = 1.0 / (double)N;

    // Outer loop over particles to maximize CPU cache locality
    for (int j = 0; j < N; j++) {
        const double w = p ? p[j] : default_weight;
        const struct reb_particle* pt = &particles[j];

        // Accumulate state vector in one contiguous pass
        m[0] += w * pt->x;
        if (M > 1) m[1] += w * pt->y;
        if (M > 2) m[2] += w * pt->z;
        if (M > 3) m[3] += w * pt->vx;
        if (M > 4) m[4] += w * pt->vy;
        if (M > 5) m[5] += w * pt->vz;
    }

    return m;
}

double** cov_particles(int N, int M, const struct reb_particle* particles, const double* p, const double* mu){
    if (N <= 0 || M <= 0 || !particles || !mu) return NULL;

    // Allocate row pointer array and continuous matrix block in two allocations
    double** S = (double**)malloc(M * sizeof(double*));
    if (!S) return NULL;

    double* block = (double*)calloc(M * M, sizeof(double));
    if (!block) {
        free(S);
        return NULL;
    }

    for (int i = 0; i < M; i++) {
        S[i] = block + (i * M);
    }

    // Default weight factor when p == NULL (unbiased sample covariance factor 1 / (N - 1))
    const bool default_p = (p == NULL);
    const double default_weight = (N > 1) ? (1.0 / (double)(N - 1)) : 1.0;

    // Accumulate outer product (dx * dx^T) into upper triangle across single pass
    for (int j = 0; j < N; j++) {
        const double w = default_p ? default_weight : p[j];
        const struct reb_particle* pt = &particles[j];

        // Extract deviation vector (x - mu) into stack memory
        double dx[6] = {0};
        dx[0] = pt->x - mu[0];
        if (M > 1) dx[1] = pt->y - mu[1];
        if (M > 2) dx[2] = pt->z - mu[2];
        if (M > 3) dx[3] = pt->vx - mu[3];
        if (M > 4) dx[4] = pt->vy - mu[4];
        if (M > 5) dx[5] = pt->vz - mu[5];

        // Symmetric outer product accumulation
        for (int r = 0; r < M; r++) {
            const double w_dx_r = w * dx[r];
            for (int c = r; c < M; c++) {
                S[r][c] += w_dx_r * dx[c];
            }
        }
    }

    // Fill symmetric lower triangle
    for (int r = 0; r < M; r++) {
        for (int c = 0; c < r; c++) {
            S[r][c] = S[c][r];
        }
    }

    return S;
}

double calc_HZ_fast_particles(int N, int dim, const struct reb_particle* particles, double* p, double* mu, double** S){
    if (N <= 0 || dim <= 0 || !particles) return 0.0;

    /* compute beta and frequently-used constants once */
    double beta = (1.0 / sqrt(2.0)) * pow(((N * (2.0 * dim + 1.0)) / 4.0), 1.0 / (dim + 4.0));
    double beta2 = beta * beta;
    double c1 = beta2 * 0.5; /* for HZ1: exp(-c1 * d) where d is Mahalanobis quadratic */
    double c2 = beta2 / (2.0 * (1.0 + beta2)); /* for HZ2 exponent */
    double two_times = 2.0 * pow((1.0 + beta2), -((double)dim) / 2.0);
    double last_term = N * pow((1.0 + 2.0 * beta2), -((double)dim) / 2.0);

    /* Compute mean and covariance directly from particle array */
    double* mu_t = mean_particles(N, dim, particles, p); /* length dim, must be freed by caller */
    double** S_t = cov_particles(N, dim, particles, p, mu_t);  /* dim x dim double** (row pointers), must be freed by caller */
    for(int i = 0; i < dim; i++){
        mu[i] = mu_t[i]; 
        for(int j = 0; j < dim; j++){
            S[i][j] = S_t[i][j]; 
        }
    }
    free(mu_t);
    if(S_t){
        free(S_t[0]); // frees the single continuous block
        free(S_t);    // frees the row pointer array
    }

    /* Cholesky of S -> R (upper triangular) */
    double* R = cholesky_upper(dim, S);
    if (!R) {
        return NAN;
    }

    /* Transform points: z_i = R^{-T} x_i (i.e. solve R^T z = x) */
    double* Z = init_v(N * dim); 
    if (!Z) {
        free(R);
        return NAN;
    }

    /* Extract particle state directly into forward substitution without staging in X */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        double temp[6]; // Stack allocation (dim <= 6) avoids per-thread heap malloc
        const struct reb_particle* pt = &particles[i];
        
        // Load particle state into local contiguous stack buffer
        temp[0] = pt->x;  temp[1] = pt->y;  temp[2] = pt->z;
        temp[3] = pt->vx; temp[4] = pt->vy; temp[5] = pt->vz;

        // Solve R^T * Z[i] = temp directly into destination memory
        solve_Rt(dim, R, temp, Z + i * dim);
    }

    /* Transform mu as well */
    double* zmu = (double*)malloc(sizeof(double) * dim);
    if (!zmu) {
        free(Z);
        free(R);
        return NAN;
    }
    solve_Rt(dim, R, mu, zmu);

    /* Compute HZ1 using symmetric pair sum */
    double HZ1 = (double)N;
    double off = 0.0;

    #pragma omp parallel for reduction(+:off) schedule(static)
    for (int i = 0; i < N; ++i) {
        const double* zi = Z + i * dim;
        for (int j = i + 1; j < N; ++j) {
            const double* zj = Z + j * dim;
            double d2 = sq_euclid(dim, zi, zj);
            off += exp(-c1 * d2);
        }
    }
    HZ1 += 2.0 * off;
    HZ1 /= (double)N;

    /* Compute HZ2: sum_i exp(-c2 * ||z_i - zmu||^2) then multiply by two_times */
    double HZ2 = 0.0;
    #pragma omp parallel for reduction(+:HZ2) schedule(static)
    for (int i = 0; i < N; ++i) {
        const double* zi = Z + i * dim;
        double d2 = sq_euclid(dim, zi, zmu);
        HZ2 += exp(-c2 * d2);
    }
    HZ2 *= two_times;

    /* Final HZ */
    double HZ = HZ1 - HZ2 + last_term;

    /* Cleanup */
    free(Z);
    free(zmu);
    free(R);

    return HZ;
}

double calc_b1_particles(int N, int dim, const struct reb_particle* particles, double* mu, double** S){
    if (N <= 0 || dim <= 0 || !particles || !mu || !S) return NAN;

    // 1. Compute mean directly into caller-provided mu
    double* mu_t = mean_particles(N, dim, particles, NULL);
    if (!mu_t) return NAN;
    for (int j = 0; j < dim; j++) {
        mu[j] = mu_t[j];
    }
    free(mu_t);

    // 2. Compute covariance matrix using cov_particles
    double** S_p = cov_particles(N, dim, particles, NULL, mu);
    if (!S_p) return NAN;

    // Copy S_p into caller S and 1D S_t buffer for LAPACKE_dposv
    double S_t[36]; // Stack allocation for dim <= 6
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            S[i][j] = S_p[i][j];
            S_t[i * dim + j] = S_p[i][j];
        }
    }

    // Free S_p (block layout: free S_p[0] then S_p)
    free(S_p[0]);
    free(S_p);

    // 3. Extract mean-centered particles directly into Xc and transpose A
    double* Xc = init_v(N * dim); 
    double* A  = init_v(dim * N); 
    if (!Xc || !A) {
        free(Xc); free(A);
        return NAN;
    }

    for (int i = 0; i < N; i++) {
        const struct reb_particle* pt = &particles[i];
        double dx[6] = {0};
        dx[0] = pt->x - mu[0];
        if (dim > 1) dx[1] = pt->y - mu[1];
        if (dim > 2) dx[2] = pt->z - mu[2];
        if (dim > 3) dx[3] = pt->vx - mu[3];
        if (dim > 4) dx[4] = pt->vy - mu[4];
        if (dim > 5) dx[5] = pt->vz - mu[5];

        for (int j = 0; j < dim; j++) {
            Xc[i * dim + j] = dx[j];
            A[j * N + i]    = dx[j];
        }
    }

    // 4. Solve S * A = Xc' via LAPACK
    int info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', dim, N, S_t, dim, A, N);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_dposv failed, info = %d\n", info);
        free(Xc); free(A);
        return NAN;
    }

    // 5. Compute Q = Xc * A (N x N matrix)
    double* Q = init_v(N * N); 
    if (!Q) {
        free(Xc); free(A);
        return NAN;
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, dim, 1.0, Xc, dim, A, N, 0.0, Q, N);

    // 6. Compute sum_ij Q_ij^3 / N^2
    double b1 = 0.0;
    #pragma omp parallel for reduction(+:b1) schedule(static)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            const double q = Q[i * N + j];
            b1 += q * q * q;
        }
    }
    b1 /= ((double)N * (double)N);

    // Cleanup
    free(Xc); free(A); free(Q);

    return b1;
}

double calc_b2_particles(int N, int dim, const struct reb_particle* particles, double* mu, double** S){
    if (N <= 0 || dim <= 0 || !particles || !mu || !S) return NAN;

    // 1. Compute mean directly into caller-provided mu
    double* mu_t = mean_particles(N, dim, particles, NULL);
    if (!mu_t) return NAN;
    for (int j = 0; j < dim; j++) {
        mu[j] = mu_t[j];
    }
    free(mu_t);

    // 2. Compute covariance matrix using cov_particles
    double** S_p = cov_particles(N, dim, particles, NULL, mu);
    if (!S_p) return NAN;

    // Copy S_p into caller S and 1D S_t buffer for LAPACKE_dposv
    double S_t[36]; // Stack allocation for dim <= 6
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            S[i][j] = S_p[i][j];
            S_t[i * dim + j] = S_p[i][j];
        }
    }

    // Free S_p (block layout: free S_p[0] then S_p)
    free(S_p[0]);
    free(S_p);

    // 3. Extract mean-centered particles directly into transpose A (dim x N)
    double* A = init_v(dim * N); 
    if (!A) return NAN;

    for (int i = 0; i < N; i++) {
        const struct reb_particle* pt = &particles[i];
        double dx[6] = {0};
        dx[0] = pt->x - mu[0];
        if (dim > 1) dx[1] = pt->y - mu[1];
        if (dim > 2) dx[2] = pt->z - mu[2];
        if (dim > 3) dx[3] = pt->vx - mu[3];
        if (dim > 4) dx[4] = pt->vy - mu[4];
        if (dim > 5) dx[5] = pt->vz - mu[5];

        for (int j = 0; j < dim; j++) {
            A[j * N + i] = dx[j];
        }
    }

    // 4. Solve S * A = Xc' via LAPACK.
    // A stores S^{-1} * Xc' as columns: A[:, i] = S^{-1} * (x_i - mu)
    int info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', dim, N, S_t, dim, A, N);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_dposv failed, info = %d\n", info);
        free(A);
        return NAN;
    }

    // 5. Compute q_ii = (x_i - mu)^T * (S^{-1} * (x_i - mu)) directly without forming full Q (N x N)
    double b2 = 0.0;
    #pragma omp parallel for reduction(+:b2) schedule(static)
    for (int i = 0; i < N; i++) {
        const struct reb_particle* pt = &particles[i];
        double dx[6] = {0};
        dx[0] = pt->x - mu[0];
        if (dim > 1) dx[1] = pt->y - mu[1];
        if (dim > 2) dx[2] = pt->z - mu[2];
        if (dim > 3) dx[3] = pt->vx - mu[3];
        if (dim > 4) dx[4] = pt->vy - mu[4];
        if (dim > 5) dx[5] = pt->vz - mu[5];

        double qii = 0.0;
        for (int j = 0; j < dim; j++) {
            qii += dx[j] * A[j * N + i]; // dot product of Xc row i and A column i
        }
        b2 += qii * qii;
    }
    b2 /= (double)N;

    // Cleanup
    free(A);

    return b2;
}

void update_moments_ukf_particles(int dim, const struct reb_particle* particles, double* W, double* mu, double** S){
    if(!particles || !W || !mu || !S || dim <= 0) return;

    int N = 2 * dim + 1;

    // 1. Compute mean across sigma points: mu = sum_i W[i] * particle_i
    for(int j = 0; j < dim; j++){
        mu[j] = 0.0;
    }

    for(int i = 0; i < N; i++){
        double w = W[i];
        const struct reb_particle* pt = &particles[i];

        mu[0] += w * pt->x;
        if (dim > 1) mu[1] += w * pt->y;
        if (dim > 2) mu[2] += w * pt->z;
        if (dim > 3) mu[3] += w * pt->vx;
        if (dim > 4) mu[4] += w * pt->vy;
        if (dim > 5) mu[5] += w * pt->vz;
    }

    // 2. Initialize S to zero
    for(int j = 0; j < dim; j++){
        for(int k = 0; k < dim; k++){
            S[j][k] = 0.0;
        }
    }

    // 3. Compute weighted covariance using weights W[i + N]
    for(int i = 0; i < N; i++){
        double w_cov = W[i + N];
        const struct reb_particle* pt = &particles[i];

        double diff[6] = {0};
        diff[0] = pt->x - mu[0];
        if (dim > 1) diff[1] = pt->y - mu[1];
        if (dim > 2) diff[2] = pt->z - mu[2];
        if (dim > 3) diff[3] = pt->vx - mu[3];
        if (dim > 4) diff[4] = pt->vy - mu[4];
        if (dim > 5) diff[5] = pt->vz - mu[5];

        for(int j = 0; j < dim; j++){
            const double w_diff_j = w_cov * diff[j];
            for (int k = j; k < dim; k++) {
                S[j][k] += w_diff_j * diff[k];
            }
        }
    }

    // Fill symmetric lower triangle
    for(int j = 0; j < dim; j++){
        for(int k = 0; k < j; k++){
            S[j][k] = S[k][j];
        }
    }

    return; 
}

double calc_NED_particles(int dim, const struct reb_particle* particles, double* W, double* mu, double** S){
    if(!particles || !W || !mu || !S || dim <= 0) return NAN;

    // 1. Update mean and covariance directly from particle subset
    update_moments_ukf_particles(dim, particles, W, mu, S);

    // 2. Compute difference between mean sigma point (particles[0]) and mu
    double diff[6] = {0};
    const struct reb_particle* p0 = &particles[0];
    diff[0] = p0->x - mu[0];
    if (dim > 1) diff[1] = p0->y - mu[1];
    if (dim > 2) diff[2] = p0->z - mu[2];
    if (dim > 3) diff[3] = p0->vx - mu[3];
    if (dim > 4) diff[4] = p0->vy - mu[4];
    if (dim > 5) diff[5] = p0->vz - mu[5];

    // 3. Flatten S into contiguous array for LAPACK Cholesky factor
    double S_t[36];
    for(int j = 0; j < dim; j++){
        for(int k = 0; k < dim; k++){
            S_t[j * dim + k] = S[j][k];
        }
    }

    // 4. Compute Cholesky decomposition S = L * L^T (Lower triangular)
    int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', dim, S_t, dim);
    if (info != 0) {
        return NAN; // Matrix is not positive-definite
    }

    // 5. Solve L * y = diff directly (forward substitution)
    //    ||y||^2 = diff^T * S^{-1} * diff = ||L^{-1} * diff||^2
    double y[6] = {0};
    for(int j = 0; j < dim; j++){
        double sum = diff[j];
        for(int k = 0; k < j; k++){
            sum -= S_t[j * dim + k] * y[k];
        }
        y[j] = sum / S_t[j * dim + j];
    }

    // 6. Return Euclidean norm of y
    double norm_sq = 0.0;
    for(int j = 0; j < dim; j++){
        norm_sq += y[j] * y[j];
    }

    return sqrt(norm_sq);
}

double* reconstruct_cut_particles(int dim, int M1, int M2, const struct reb_particle* particles, const double* W){
    if (dim <= 0 || !particles || !W) return NULL;

    if (M1 < 1 || M1 > M2) {
        fprintf(stderr, "M1 must be an integer between 1 and the CUT order M2.\n");
        return NULL;
    }

    // Determine number of CUT sigma points N
    int N = 0;
    if (M2 == 4) {
        N = 2 * dim + (1 << dim) + 1;
    } else if (M2 == 6) {
        if (dim <= 6) {
            N = 2 * dim * dim + (1 << dim) + 1;
        } else {
            N = 2 * dim + (1 << dim) + (4 * dim * (dim - 1) * (dim - 2) / 3) + 1;
        }
    } else {
        fprintf(stderr, "Unsupported CUT order M2 = %d\n", M2);
        return NULL;
    }

    // Compute size = dim^M1
    int size = 1;
    for (int i = 0; i < M1; i++) size *= dim;

    double* M = (double*)calloc(size, sizeof(double));
    if (!M) return NULL;

    // Standardize/normalize weight sum into stack buffer or flat vector
    double wsum = 0.0;
    for (int i = 0; i < N; i++) wsum += W[i];

    // Stack allocation for weights if N <= 256, else heap
    double W_norm_stack[256];
    double* W_norm = (N <= 256) ? W_norm_stack : (double*)malloc(N * sizeof(double));
    if (!W_norm) { free(M); return NULL; }

    for (int i = 0; i < N; i++) {
        W_norm[i] = W[i] / wsum;
    }

    // 1. Compute Mean Vector mu (Length dim)
    double mu[6] = {0}; // Stack allocation for dim <= 6
    for (int n = 0; n < N; n++) {
        const double w = W_norm[n];
        const struct reb_particle* pt = &particles[n];

        mu[0] += w * pt->x;
        if (dim > 1) mu[1] += w * pt->y;
        if (dim > 2) mu[2] += w * pt->z;
        if (dim > 3) mu[3] += w * pt->vx;
        if (dim > 4) mu[4] += w * pt->vy;
        if (dim > 5) mu[5] += w * pt->vz;
    }

    // If M1 == 1, return mean vector directly
    if (M1 == 1) {
        for (int i = 0; i < dim; i++) {
            M[i] = mu[i];
        }
        if (W_norm != W_norm_stack) free(W_norm);
        return M;
    }

    // 2. Case M1 == 2: Covariance Matrix (dim x dim, flattened)
    if (M1 == 2) {
        for (int n = 0; n < N; n++) {
            const double w = W_norm[n];
            const struct reb_particle* pt = &particles[n];

            double dx[6] = {0};
            dx[0] = pt->x - mu[0];
            if (dim > 1) dx[1] = pt->y - mu[1];
            if (dim > 2) dx[2] = pt->z - mu[2];
            if (dim > 3) dx[3] = pt->vx - mu[3];
            if (dim > 4) dx[4] = pt->vy - mu[4];
            if (dim > 5) dx[5] = pt->vz - mu[5];

            for (int i = 0; i < dim; i++) {
                const double w_dxi = w * dx[i];
                for (int j = i; j < dim; j++) {
                    M[i * dim + j] += w_dxi * dx[j];
                }
            }
        }

        // Mirror upper triangle to lower triangle
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < i; j++) {
                M[i * dim + j] = M[j * dim + i];
            }
        }

        if (W_norm != W_norm_stack) free(W_norm);
        return M;
    }

    // 3. Case M1 > 2: High-Order Outer Product Accumulation
    for (int n = 0; n < N; n++) {
        const double w = W_norm[n];
        const struct reb_particle* pt = &particles[n];

        double dx[6] = {0};
        dx[0] = pt->x - mu[0];
        if (dim > 1) dx[1] = pt->y - mu[1];
        if (dim > 2) dx[2] = pt->z - mu[2];
        if (dim > 3) dx[3] = pt->vx - mu[3];
        if (dim > 4) dx[4] = pt->vy - mu[4];
        if (dim > 5) dx[5] = pt->vz - mu[5];

        // Compute outer_M on stack buffer to eliminate per-particle malloc
        double* T = outer_M(dx, dim, M1);
        if (T) {
            for (int j = 0; j < size; j++) {
                M[j] += w * T[j];
            }
            free(T);
        }
    }

    if (W_norm != W_norm_stack) free(W_norm);

    return M;
}

double calc_bt1_particles(int N, int dim, int M, const struct reb_particle* particles, double* W, double* mu, double** S){
    if (N <= 0 || dim <= 0 || !particles || !W || !mu || !S) return NAN;

    // 1. Reconstruct mean using particle-native CUT function
    double* mu_t = reconstruct_cut_particles(dim, 1, M, particles, W);
    if (!mu_t) return NAN;
    for (int j = 0; j < dim; j++) {
        mu[j] = mu_t[j];
    }
    free(mu_t);

    // 2. Reconstruct covariance using particle-native CUT function
    double* S_t = reconstruct_cut_particles(dim, 2, M, particles, W);
    if (!S_t) return NAN;

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            S[i][j] = S_t[i * dim + j];
        }
    }

    // 3. Extract mean-centered particles directly into Zc and transposed X
    double* Zc = (double*)malloc(sizeof(double) * N * dim);
    double* X  = (double*)malloc(sizeof(double) * dim * N);
    if (!Zc || !X) {
        free(S_t); free(Zc); free(X);
        return NAN;
    }

    for (int i = 0; i < N; i++) {
        const struct reb_particle* pt = &particles[i];
        double dx[6] = {0};
        dx[0] = pt->x - mu[0];
        if (dim > 1) dx[1] = pt->y - mu[1];
        if (dim > 2) dx[2] = pt->z - mu[2];
        if (dim > 3) dx[3] = pt->vx - mu[3];
        if (dim > 4) dx[4] = pt->vy - mu[4];
        if (dim > 5) dx[5] = pt->vz - mu[5];

        for (int j = 0; j < dim; j++) {
            Zc[i * dim + j] = dx[j];
            X[j * N + i]    = dx[j];
        }
    }

    // 4. Solve S * X = Zc' via LAPACK (dposv overwrites S_t and X)
    int info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', dim, N, S_t, dim, X, N);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_dposv failed, info = %d\n", info);
        free(S_t); free(Zc); free(X);
        return NAN;
    }

    // 5. Compute Q = Zc * X (N x N matrix)
    double* Q = (double*)malloc(sizeof(double) * N * N);
    if (!Q) {
        free(S_t); free(Zc); free(X);
        return NAN;
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, dim, 1.0, Zc, dim, X, N, 0.0, Q, N);

    // 6. Compute bt1 = sum_ij W_i * W_j * Q_ij^3
    double bt1 = 0.0;
    #pragma omp parallel for reduction(+:bt1) schedule(static)
    for (int i = 0; i < N; i++) {
        double wi = W[i];
        for (int j = 0; j < N; j++) {
            double wij = wi * W[j];
            const double q = Q[i * N + j];
            bt1 += wij * q * q * q;
        }
    }

    // Cleanup
    free(S_t); free(Zc); free(X); free(Q);

    return bt1;
}

double calc_bt2_particles(int N, int dim, int M, const struct reb_particle* particles, double* W, double* mu, double** S){
    if (N <= 0 || dim <= 0 || !particles || !W || !mu || !S) return NAN;

    // 1. Reconstruct mean using particle-native CUT function
    double* mu_t = reconstruct_cut_particles(dim, 1, M, particles, W);
    if (!mu_t) return NAN;
    for (int j = 0; j < dim; j++) {
        mu[j] = mu_t[j];
    }
    free(mu_t);

    // 2. Reconstruct covariance using particle-native CUT function
    double* S_t = reconstruct_cut_particles(dim, 2, M, particles, W);
    if (!S_t) return NAN;

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            S[i][j] = S_t[i * dim + j];
        }
    }

    // 3. Extract mean-centered particles directly into transposed X (dim x N)
    double* X = (double*)malloc(sizeof(double) * dim * N);
    if (!X) {
        free(S_t);
        return NAN;
    }

    for (int i = 0; i < N; i++) {
        const struct reb_particle* pt = &particles[i];
        double dx[6] = {0};
        dx[0] = pt->x - mu[0];
        if (dim > 1) dx[1] = pt->y - mu[1];
        if (dim > 2) dx[2] = pt->z - mu[2];
        if (dim > 3) dx[3] = pt->vx - mu[3];
        if (dim > 4) dx[4] = pt->vy - mu[4];
        if (dim > 5) dx[5] = pt->vz - mu[5];

        for (int j = 0; j < dim; j++) {
            X[j * N + i] = dx[j];
        }
    }

    // 4. Solve S * X = Zc' via LAPACK.
    // X now holds S^{-1} * Zc' as columns: X[:, i] = S^{-1} * (z_i - mu)
    int info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', dim, N, S_t, dim, X, N);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_dposv failed, info = %d\n", info);
        free(S_t); free(X);
        return NAN;
    }

    // 5. Compute q_ii = (z_i - mu)^T * S^{-1} * (z_i - mu) directly
    //    Eliminates forming full Q (N x N) matrix via GEMM.
    double bt2 = 0.0;
    #pragma omp parallel for reduction(+:bt2) schedule(static)
    for (int i = 0; i < N; i++) {
        const struct reb_particle* pt = &particles[i];
        double dx[6] = {0};
        dx[0] = pt->x - mu[0];
        if (dim > 1) dx[1] = pt->y - mu[1];
        if (dim > 2) dx[2] = pt->z - mu[2];
        if (dim > 3) dx[3] = pt->vx - mu[3];
        if (dim > 4) dx[4] = pt->vy - mu[4];
        if (dim > 5) dx[5] = pt->vz - mu[5];

        double qii = 0.0;
        for (int j = 0; j < dim; j++) {
            qii += dx[j] * X[j * N + i]; // Dot product of row i (Zc) and column i (X)
        }
        bt2 += W[i] * qii * qii;
    }

    // Cleanup
    free(S_t); free(X);

    return bt2;
}