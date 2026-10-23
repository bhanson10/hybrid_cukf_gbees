// cukf.c
// Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

#include "cukf.h"

double** eye(int dim){

    double** Id = init_m(dim, dim); 
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            if(i == j){
                Id[i][j] = 1.0;
            }else{
                Id[i][j] = 0.0;
            }
        }
    }

    return Id; 
}

static void pack_matrix(double** A, double* out, int n){
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            out[i*n + j] = A[i][j];
}

static void unpack_matrix(double* A, double** out, int n){
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            out[i][j] = A[i*n + j];
}

static void matmul(int n, const double* A, const double* B, double* C){
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            C[i*n + j] = 0.0;
            for (int k = 0; k < n; k++)
                C[i*n + j] += A[i*n + k] * B[k*n + j];
        }
}

static void transpose2(int n, const double* A, double* AT){
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            AT[j*n + i] = A[i*n + j];
}

static double sylv_scalar(double a, double b, double c){
    return c / (a + b);
}


static int sqrtm_quasitri(int n, const double* T, double* R){
    memset(R, 0, sizeof(double)*n*n);

    for (int i = 0; i < n; i++)
    {
        if (T[i*n + i] < 0.0)
            return -1;

        R[i*n + i] = sqrt(T[i*n + i]);

        for (int j = i+1; j < n; j++)
        {
            double s = T[i*n + j];

            for (int k = i+1; k < j; k++)
                s -= R[i*n + k] * R[k*n + j];

            R[i*n + j] =
                sylv_scalar(R[i*n + i], R[j*n + j], s);
        }
    }

    return 0;
}

double** sqrtm(int dim, double** P){
    int info;

    /* Allocate result */
    double** X = init_m(dim, dim);

    /* Contiguous working arrays */
    double* A  = init_v(dim * dim);
    double* Q  = init_v(dim * dim);
    double* T  = init_v(dim * dim);
    double* R  = init_v(dim * dim);
    double* QT = init_v(dim * dim);
    double* tmp= init_v(dim * dim);

    double* wr = init_v(dim);
    double* wi = init_v(dim);
    lapack_int sdim;

    /* Pack input */
    pack_matrix(P, A, dim);

    /* ----- Real Schur decomposition ----- */
    info = LAPACKE_dgees(
        LAPACK_ROW_MAJOR,
        'V',
        'N',
        NULL,
        dim,
        A,
        dim,
        &sdim,
        wr,
        wi,
        Q,
        dim
    );

    if (info != 0)
        return NULL;

    memcpy(T, A, sizeof(double)*dim*dim);

    /* ----- sqrt(T) ----- */
    if (sqrtm_quasitri(dim, T, R) != 0)
        return NULL;

    /* ----- Reconstruct sqrt(A) ----- */
    transpose2(dim, Q, QT);

    matmul(dim, Q, R, tmp);
    matmul(dim, tmp, QT, A);

    /* Unpack result */
    unpack_matrix(A, X, dim);

    /* Free */
    free(A); free(Q); free(T); free(R);
    free(QT); free(tmp); free(wr); free(wi);

    return X;
}

void cut4(int dim, double* Z, double* W, double* mu, double** S){

    const double r[2][2] = {{1.4861736616297834, 2.6060099476935847}, 
                            {3.2530871022700643, 1.190556300661233}}; 
    const double w[3][2] = {{0.5811010092660772, 0.41553535186548973}, 
                            {0.20498484723245053, 0.021681819434216532}, 
                            {0.00446464813451093, 0.12443434259941118}};

    int N = 2 * dim + (1 << dim) + 1; 
    double** s = principal_points(dim); 
    double** cN = conj_N_points(dim);

    // Zero point
    for (int j = 0; j < dim; j++) Z[0 * dim + j] = 0.0;
    W[0] = w[0][dim - 1];
    int last_idx = 0;

    // Principal axis points
    for (int i = 0; i < 2 * dim; i++){
        for (int j = 0; j < dim; j++){
            Z[last_idx + 1 + i * dim + j] = r[0][dim - 1] * s[i][j];
        }
        W[last_idx + 1 + i] = w[1][dim - 1];
    }
    last_idx += 2 * dim;

    // Conjugate axis points
    int nCN = 1 << dim;
    for (int i = 0; i < nCN; i++) {
        for (int j = 0; j < dim; j++) {
            Z[last_idx + 1 + i * dim + j] = r[1][dim - 1] * cN[i][j];
        }
        W[last_idx + 1 + i] = w[2][dim - 1];
    }

    double** L = sqrtm(dim, S);  // sqrtm(S)

    for (int i = 0; i < N; i++) {
        double* tmp = mvtimes(dim, dim, L, &Z[i * dim]); 
        for(int j = 0; j < dim; j++){
            Z[i * dim + j] = tmp[j] + mu[j]; 
        }
        free(tmp);
    }

    // Free intermediate arrays
    for (int i = 0; i < dim; i++) free(L[i]);
    free(L);
    for (int i = 0; i < dim; i++) free(s[i]);
    free(s);
    for (int i = 0; i < (1 << dim); i++) free(cN[i]);
    free(cN);

    return; 
}

void cut6(int dim, double* Z, double* W, double* mu, double** S){
    // Constants
    const double r[3][8] = {
        {2.4494897427, 2.3587090379, 2.2520650012, 2.1213203430, 1.9488352799, 2.5512003554, 2.4494897427, 2.3439068809},
        {1.1147379454, 1.1198362859, 1.1260325006, 1.1338934189, 1.1445968942, 0.9642630978, 0.9999999999, 1.0232621870},
        {3.2004125801, 3.1421303838, 3.0763780026, 3.0000000003, 2.9068006056, 2.3255766977, 2.4494897428, 2.5342865662}
    };
    const double w[4][8] = {
        {0.3658770242, 0.3124789714, 0.2499999999, 0.1728395055, 0.0674637121, 0.0896488470, 0.0740740740, 0.0421898530},
        {0.0277777777, 0.0290351301, 0.0306601632, 0.0329218107, 0.0365072564, 0.0126940628, 0.0138888888, 0.0150764081},
        {0.1302876649, 0.0633844605, 0.0306601632, 0.0147033607, 0.0069487173, 0.0048594459, 0.0023437500, 0.0011342720},
        {0.0004653012, 0.0005195469, 0.0005898367, 0.0006858710, 0.0008288549, 0.0003950899, 0.0002314814, 0.0001572730}
    };

    int N;
    if(dim <= 6){
        N = 2 * dim * dim + (1 << dim) + 1;
    }else{
        N = 2 * dim + (1 << dim) + (4 * dim * (dim - 1) * (dim - 2) / 3) + 1;
    }

    // Generate principal and conjugate points
    double** s = principal_points(dim);   // 2*dim x dim
    double** cN = conj_N_points(dim);     // 2^dim x dim

    // Zero point
    for(int j = 0; j < dim; j++) Z[0 * dim + j] = 0.0;
    W[0] = w[0][dim - 2];
    int last_idx = 0;

    // Principal axis points
    for (int i = 0; i < 2 * dim; i++){
        for (int j = 0; j < dim; j++){
            Z[(last_idx + 1 + i) * dim + j] = r[0][dim - 2] * s[i][j];
        }
        W[last_idx + 1 + i] = w[1][dim - 2];
    }
    last_idx += 2 * dim;

    // Conjugate N-axis points
    int nCN = 1 << dim;
    for(int i = 0; i < nCN; i++) {
        for(int j = 0; j < dim; j++)
            Z[(last_idx + i + 1) * dim + j] = r[1][dim - 2] * cN[i][j];
        W[last_idx + i + 1] = w[2][dim - 2];
    }
    last_idx += nCN;

    // Conjugate M-axis points
    double** cM = NULL;
    int mCN;
    if(dim <= 6) {
        cM = conj_M_points(dim, 2);
        mCN = 2 * dim * (dim - 1);
    } else {
        cM = conj_M_points(dim, 3);
        mCN = 4 * dim * (dim - 1) * (dim - 2) / 3;
    }

    for(int i = 0; i < mCN; i++) {
        for(int j = 0; j < dim; j++)
            Z[(last_idx + i + 1) * dim + j] = r[2][dim - 2] * cM[i][j];
        W[last_idx + i + 1] = w[3][dim - 2];
    }
    last_idx += mCN;

    // Free intermediate arrays
    for(int i = 0; i < dim; i++) free(s[i]);
    free(s);
    for(int i = 0; i < (1 << dim); i++) free(cN[i]);
    free(cN);
    for(int i = 0; i < mCN; i++) free(cM[i]);
    free(cM);
    
    double** L = sqrtm(dim, S); // sqrtm(S)

    for (int i = 0; i < N; i++) {
        double* tmp = mvtimes(dim, dim, L, &Z[i * dim]); 
        for(int j = 0; j < dim; j++){
            Z[i * dim + j] = tmp[j] + mu[j]; 
        }
        free(tmp);
    }

    for(int i = 0; i < dim; i++) free(L[i]);
    free(L);
}

double** principal_points(int dim){

    double** X = init_m(2 * dim, dim); 
    double** Id = eye(dim); 
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            X[i][j] = Id[i][j];
        }
    }
    for(int i = dim; i < 2 * dim; i++){
        for(int j = 0; j < dim; j++){
            X[i][j] = -Id[i - dim][j]; 
        }
    }

    return X; 
}

double** conj_N_points(int dim){

    if (dim == 1) {
        double** X = init_m(2, 1);
        X[0][0] = 1.0;
        X[1][0] = -1.0;
        return X;
    }

    int half = 1 << (dim - 1);  // 2^(dim-1)
    int full = 1 << dim;        // 2^dim

    // recursive call
    double** D = conj_N_points(dim - 1);  // D is 2^(dim-1) x (dim-1)

    double** X = init_m(full, dim);

    // Top half: first column = 1, rest from D
    for (int i = 0; i < half; i++) {
        X[i][0] = 1.0;
        for (int j = 1; j < dim; j++) {
            X[i][j] = D[i][j - 1];
        }
    }

    // Bottom half: negate top half
    for (int i = 0; i < half; i++) {
        for (int j = 0; j < dim; j++) {
            X[i + half][j] = -X[i][j];
        }
    }

    // Free D
    for (int i = 0; i < half; i++) {
        free(D[i]);
    }
    free(D);

    return X;
}

double** conj_M_points(int dim, int order){
    // Generate conjugate M-axis points for given dim and order
    // Output: 2^order * (dim choose order) x dim

    double** D = conj_N_points(order);  // 2^order x order

    int comb_count = nchoosek(dim, order);  // number of combinations
    int pow2M = 1 << order;                 // 2^order

    // Allocate output and zero-initialize
    double** X = init_m(pow2M * comb_count, dim);
    for(int i = 0; i < pow2M * comb_count; i++)
        for(int j = 0; j < dim; j++)
            X[i][j] = 0.0;

    // Generate first combination: 0,1,...,order-1
    int* indices = (int*)malloc(order * sizeof(int));
    for(int i = 0; i < order; i++) indices[i] = i;

    int cnt = 0;  // row offset in X
    while(1) {
        // Copy D into X at positions specified by indices
        for(int row = 0; row < pow2M; row++) {
            for(int j = 0; j < order; j++) {
                X[cnt + row][indices[j]] = D[row][j];
            }
        }
        cnt += pow2M;

        // Generate next combination in lexicographic order
        int i;
        for(i = order - 1; i >= 0; i--){
            if(indices[i] < dim - order + i){
                indices[i]++;
                for(int j = i + 1; j < order; j++)
                    indices[j] = indices[j - 1] + 1;
                break;
            }
        }
        if(i < 0) break;  // finished all combinations
    }

    // Free temporary arrays
    int D_rows = 1 << order;
    for(int i = 0; i < D_rows; i++) free(D[i]);
    free(D);
    free(indices);

    return X;
}

int nchoosek(int N, int M){
    if (M > N) return 0;
    if (M == 0 || M == N) return 1;
    int res = 1;
    for (int i = 1; i <= M; i++) {
        res = res * (N - i + 1) / i;
    }
    return res;
}

double* reconstruct_cut(int dim, int M1, int M2, double* Z, double* W){

    int size = 1;
    for (int i = 0; i < M1; i++) size *= dim;   
    double* M = init_v(size); 

    int N; 
    if(M2 == 4){
        N = 2 * dim + (1 << dim) + 1; 
    }else if(M2 == 6){
        if(dim <= 6){
            N = 2 * dim * dim + (1 << dim) + 1;
        }else{
            N = 2 * dim + (1 << dim) + (4 * dim * (dim - 1) * (dim - 2) / 3) + 1;
        }
    }

    if((M1 < 1)||(M1 > M2)){
        fprintf(stderr, "M must be an integer between 1 and the CUT order.");
        return NULL; 
    }

    double wsum = 0.0;
    for (int i = 0; i < N; i++) wsum += W[i];
    for (int i = 0; i < N; i++) W[i] /= wsum;

    double* mu = init_v(dim); 
    for(int i = 0; i < dim; i++){
        mu[i] = 0;
        for(int n = 0; n < N; n++){
            mu[i] += Z[n * dim + i] * W[n];
        }
    }

    if(M1 == 1){
        for(int i = 0; i < dim; i++){
            M[i] = mu[i]; 
        }

        return M; 
    }

    double** Zc = init_m(N, dim); 
    double** Zcw = init_m(N, dim); 
    for(int i = 0; i < dim; i++){
        for(int n = 0; n < N; n++){
            Zc[n][i]  = Z[n * dim + i] - mu[i]; 
            Zcw[n][i] = Zc[n][i] * W[n]; 
        }
    }

    if(M1 == 2){
        double** ZcwT = transpose(N, dim, Zcw);
        double** S = mmtimes(dim, N, dim, ZcwT, Zc); 
        int count = 0; 
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                M[count++] = S[i][j]; 
            }
        }

        for(int i = 0; i < N; i++){
            free(Zc[i]);
            free(Zcw[i]);
        }
        free(Zc); free(Zcw);
        for(int i = 0; i < dim; i++){
            free(ZcwT[i]); 
            free(S[i]); 
        }
        free(ZcwT); free(S); 

        return M; 
    }

    double* v = NULL; 
    double* T = NULL; 
    for(int i = 0; i < N; i++){
        v = Zc[i]; 
        T = outer_M(v, dim, M1); 
        for(int j = 0; j < size; j++){
            M[j] += W[i] * T[j]; 
        }
    }
    free(v); free(T); 

    return M; 
}

double* outer_M(double* v, int d, int M){

    if(M == 1){
        double* T = init_v(d); 
        for(int i=0; i<d; i++) T[i] = v[i];
        return T;
    }else{

        double* Tprev = outer_M(v, d, M - 1); 
        int prev_size = 1;
        for (int i = 0; i < M - 1; i++) prev_size *= d;

        int curr_size = prev_size * d;
        double* T = init_v(curr_size);

        for (int i = 0; i < prev_size; i++) {
            for (int j = 0; j < d; j++) {
                T[i * d + j] = Tprev[i] * v[j];
            }
        }

        free(Tprev);
        return T;
    }
}

void measurement_update_cukf(void (*h)(double t, double* x, double* y, double* coef), double** Z, double** W, double* y, double* coef, double** Sy, double t, int N, int dim, int y_dim, int M){
    
    // ------------------------------------------------------------
    // Reconstruct state moments
    // ------------------------------------------------------------
    double* mu = reconstruct_cut(dim, 1, M, *Z, *W);
    double* S_t = reconstruct_cut(dim, 2, M, *Z, *W);
    double** S = init_m(dim, dim);
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            S[i][j] = S_t[i * dim + j];
        }
    }
    free(S_t);

    // ------------------------------------------------------------
    // Transform CUT points through measurement model
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
        y_mu[j] = 0.0;
        for(int i = 0; i < N; i++){
            y_mu[j] += (*W)[i] * Y[i][j];
        }
    }

    // ------------------------------------------------------------
    // Measurement covariance
    // ------------------------------------------------------------

    double** P_yy = init_m(y_dim, y_dim);
    for(int i = 0; i < y_dim; i++){
        for(int j = 0; j < y_dim; j++){
            P_yy[i][j] = 0.0;
        }
    }
    double** P_xy = init_m(dim, y_dim);
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
                P_yy[j][l] += (*W)[i] * dy[j] * dy[l];
            }
        }

        // P_xy
        for(int j = 0; j < y_dim; j++){
            for(int l = 0; l < dim; l++){
                double dx = (*Z)[i * dim + l] - mu[l];
                P_xy[l][j] += (*W)[i] * dx * dy[j];
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
    // State update
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
    // Regenerate CUT
    // ------------------------------------------------------------
    if(M == 4){
        cut4(dim, *Z, *W, mu, S);
    }else if(M == 6){
        cut6(dim, *Z, *W, mu, S);
    }

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


void measurement_update_cukf_pme(void (*h)(double t, double* x, double* y, double* coef), double** Z, double** W, double* X, double* P, int num_pts, double* y, double* coef, double** Sy, double t, int N, int dim, int y_dim, int M){
    
    // printf("\nWriting prior...\n");
    // write_pdf("../../../results/PCR3BP/testing/prior.txt", num_pts, dim, X, P);

    // ------------------------------------------------------------
    // Calculate the predicted measurement for each point and its probability
    // ------------------------------------------------------------
    double* xi = init_v(dim);
    double* yi = init_v(y_dim);
    double Pi; 
    double sum = 0; 
    for(int i = 0; i < num_pts; i++){
        for(int j = 0; j < dim; j++){
            xi[j] = X[i * dim + j]; 
        }
        h(t, xi, yi, coef);
        Pi = gauss_probability(y_dim, yi, y, Sy);
        P[i] *= Pi; 
        sum += P[i];
    }
    for(int i = 0; i < num_pts; i++) P[i] /= sum; 

    // printf("\nWriting posterior...\n");
    // write_pdf("../../../results/PCR3BP/testing/posterior.txt", num_pts, dim, X, P);
    
    // fit gaussian to posterior grid
    double* mu = init_v(dim);
    double** S = init_m(dim, dim);
    weighted_mean_cov_2(X, P, num_pts, dim, mu, S);

    // ------------------------------------------------------------
    // Regenerate CUT
    // ------------------------------------------------------------
    if(M == 4){
        cut4(dim, *Z, *W, mu, S);
    }else if(M == 6){
        cut6(dim, *Z, *W, mu, S);
    }

    free(xi); 
    free(yi); 

    return;
}

void calc_bt1_bt2(int N, int dim, int M, double* Z, double* W, double* mu, double** S, double* bt1, double* bt2){

    double* mu_t = reconstruct_cut(dim, 1, M, Z, W);
    double* S_t  = reconstruct_cut(dim, 2, M, Z, W);
    for(int i = 0; i < dim; i++){
        mu[i] = mu_t[i]; 
        for(int j = 0; j < dim; j++){
            S[i][j] = S_t[i * dim + j]; 
        }
    }

    double* Zc = init_v(N * dim);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            Zc[i * dim + j] = Z[i * dim + j] - mu_t[j];
        }
    }


    double* X = init_v(dim * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            X[j * N + i] = Zc[i * dim + j];
        }
    }

    /*
     * -------------------------------------------------------
     *
     * dposv overwrites A and B.
     *
     * S * X = Zc'
     * 
     * -------------------------------------------------------
     */
    int info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', dim, N, S_t, dim, X, N);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_dposv failed, info = %d\n", info);
        free(mu_t);
        free(S_t);
        free(Zc);
        free(X);
        *bt1 = NAN;
        *bt2 = NAN;
        return;
    }

    /*
     * -------------------------------------------------------
     *
     * Q = Zc * Sinv_ZcT;
     *
     *     (N x dim) * (dim x N)
     *
     *     = N x N
     * -------------------------------------------------------
     */
    double* Q =init_v(N * N); 
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, dim, 1.0, Zc, dim, X, N, 0.0, Q, N);

    /*
     * -------------------------------------------------------
     *
     * q_diag = diag(Q);
     * -------------------------------------------------------
     */

    /*
     * -------------------------------------------------------
     *
     * W_outer = W * W';
     *
     * bt1 = sum(sum(W_outer .* (Q.^3)));
     *
     * Rather than explicitly forming W_outer, calculate
     *
     * bt1 = sum_ij W_i W_j Q_ij^3
     *
     * directly.
     * -------------------------------------------------------
     */

    *bt1 = 0.0;
    for (int i = 0; i < N; i++) {
        const double wi = W[i];
        for (int j = 0; j < N; j++) {
            const double wij = wi * W[j];
            const double q = Q[i*N + j];
            *bt1 += wij * q * q * q;
        }
    }

    /*
     * -------------------------------------------------------
     * MATLAB:
     *
     * bt2 = sum(W .* (q_diag.^2));
     * -------------------------------------------------------
     */

    *bt2 = 0.0;
    for (int i = 0; i < N; i++){
        const double qii = Q[i * N + i];
        *bt2 += W[i] * qii * qii;
    }

    free(mu_t);
    free(S_t);
    free(Zc);
    free(X);
    free(Q);

    return;
}


double calc_bt1(int N, int dim, int M, double* Z, double* W, double* mu, double** S){

    double* mu_t = reconstruct_cut(dim, 1, M, Z, W);
    double* S_t  = reconstruct_cut(dim, 2, M, Z, W);
    for(int i = 0; i < dim; i++){
        mu[i] = mu_t[i]; 
        for(int j = 0; j < dim; j++){
            S[i][j] = S_t[i * dim + j]; 
        }
    }

    double* Zc = init_v(N * dim);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            Zc[i * dim + j] = Z[i * dim + j] - mu_t[j];
        }
    }


    double* X = init_v(dim * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            X[j * N + i] = Zc[i * dim + j];
        }
    }

    /*
     * -------------------------------------------------------
     *
     * dposv overwrites A and B.
     *
     * S * X = Zc'
     * 
     * -------------------------------------------------------
     */
    int info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', dim, N, S_t, dim, X, N);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_dposv failed, info = %d\n", info);
        free(mu_t);
        free(S_t);
        free(Zc);
        free(X);
        return NAN;
    }

    /*
     * -------------------------------------------------------
     *
     * Q = Zc * Sinv_ZcT;
     *
     *     (N x dim) * (dim x N)
     *
     *     = N x N
     * -------------------------------------------------------
     */
    double* Q =init_v(N * N); 
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, dim, 1.0, Zc, dim, X, N, 0.0, Q, N);

    /*
     * -------------------------------------------------------
     *
     * q_diag = diag(Q);
     * -------------------------------------------------------
     */

    /*
     * -------------------------------------------------------
     *
     * W_outer = W * W';
     *
     * bt1 = sum(sum(W_outer .* (Q.^3)));
     *
     * Rather than explicitly forming W_outer, calculate
     *
     * bt1 = sum_ij W_i W_j Q_ij^3
     *
     * directly.
     * -------------------------------------------------------
     */

    double bt1 = 0.0;
    for (int i = 0; i < N; i++) {
        const double wi = W[i];
        for (int j = 0; j < N; j++) {
            const double wij = wi * W[j];
            const double q = Q[i*N + j];
            bt1 += wij * q * q * q;
        }
    }

    free(mu_t);
    free(S_t);
    free(Zc);
    free(X);
    free(Q);

    return bt1;
}


double calc_bt2(int N, int dim, int M, double* Z, double* W, double* mu, double** S){

    double* mu_t = reconstruct_cut(dim, 1, M, Z, W);
    double* S_t  = reconstruct_cut(dim, 2, M, Z, W);
    for(int i = 0; i < dim; i++){
        mu[i] = mu_t[i]; 
        for(int j = 0; j < dim; j++){
            S[i][j] = S_t[i * dim + j]; 
        }
    }

    double* Zc = init_v(N * dim);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            Zc[i * dim + j] = Z[i * dim + j] - mu_t[j];
        }
    }


    double* X = init_v(dim * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            X[j * N + i] = Zc[i * dim + j];
        }
    }

    /*
     * -------------------------------------------------------
     *
     * dposv overwrites A and B.
     *
     * S * X = Zc'
     * 
     * -------------------------------------------------------
     */
    int info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', dim, N, S_t, dim, X, N);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_dposv failed, info = %d\n", info);
        free(mu_t);
        free(S_t);
        free(Zc);
        free(X);
        return NAN;
    }

    /*
     * -------------------------------------------------------
     *
     * Q = Zc * Sinv_ZcT;
     *
     *     (N x dim) * (dim x N)
     *
     *     = N x N
     * -------------------------------------------------------
     */
    double* Q =init_v(N * N); 
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, dim, 1.0, Zc, dim, X, N, 0.0, Q, N);

    /*
     * -------------------------------------------------------
     * MATLAB:
     *
     * bt2 = sum(W .* (q_diag.^2));
     * -------------------------------------------------------
     */

    double bt2 = 0.0;
    for (int i = 0; i < N; i++){
        const double qii = Q[i * N + i];
        bt2 += W[i] * qii * qii;
    }

    free(mu_t);
    free(S_t);
    free(Zc);
    free(X);
    free(Q);

    return bt2;
}