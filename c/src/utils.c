// utils.c
// Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

#include "utils.h"

double* init_v(int N){
    double* v = (double*)malloc(N * sizeof(double));
    if (v == NULL) {
        perror("init_v: row malloc failed");
        return NULL;
    }

    return v; 
}

int* init_v_i(int N){
    int* v = (int*)malloc(N * sizeof(int));
    if (v == NULL) {
        perror("init_v: row malloc failed");
        return NULL;
    }

    return v; 
}

double** init_m(int N, int M){
    double** m = (double**)malloc(N * sizeof(double*));
    if (m == NULL) {
        perror("init_m: row malloc failed");
        return NULL;
    }
    for (int i = 0; i < N; i++) {
        m[i] = (double*)malloc(M * sizeof(double));
        if (m[i] == NULL) {
            perror("init_m: column malloc failed");
            for (int k = 0; k < i; k++) free(m[k]);
            free(m);
            return NULL;
        }
    }

    return m; 
}

int** init_m_i(int N, int M){
    int** m = (int**)malloc(N * sizeof(int*));
    if (m == NULL) {
        perror("init_m: row malloc failed");
        return NULL;
    }
    for (int i = 0; i < N; i++) {
        m[i] = (int*)malloc(M * sizeof(int));
        if (m[i] == NULL) {
            perror("init_m: column malloc failed");
            for (int k = 0; k < i; k++) free(m[k]);
            free(m);
            return NULL;
        }
    }

    return m; 
}

double** cast_matrix(int N, int M, double matrix[N][M]){

    double** result = init_m(N, M); 
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            result[i][j] = matrix[i][j];
        }
    }

    return result; 
}

double* mvtimes(int N, int M, double** matrix_a, const double* vector_b){

    double* result = init_v(N); 

    for (int i = 0; i < N; i++) {
        result[i] = 0.0;
        for (int j = 0; j < M; j++) {
            result[i] += matrix_a[i][j] * vector_b[j];
        }
    }
    return result;
}

double** mmtimes(int N, int M, int L, double** matrix_a, double** matrix_b){
    double** result = init_m(N, L);
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < L; j++) {
            result[i][j] = 0;
            for (int k = 0; k < M; k++) {
                result[i][j] += matrix_a[i][k] * matrix_b[k][j];
            }
        }
    }

    return result; 
}

void print_vector(int N, double* vector){
    printf("\n["); 
    for(int i = 0; i < N - 1; i++){
        printf("%.15e ", vector[i]); 
    }
    printf("%.15e]\n", vector[N - 1]); 

    return; 
}

void print_vector_T(int N, double* vector){
    printf("\n"); 
    for(int i = 0; i < N; i++){
        printf("%.15e\n", vector[i]); 
    }   
    printf("\n"); 
    return; 
}

void print_matrix(int N, int M, double** matrix){
    printf("\n"); 
    for(int i = 0; i < N; i++){
        if(i == 0){
            printf("[");
        }else{
            printf(" ");
        }
        for(int j = 0; j < M - 1; j++){
            if(matrix[i][j] < 0){
                printf("%.15e ", matrix[i][j]); 
            }else{
                printf("%.16e ", matrix[i][j]); 
            }
            
        }
        if(i == (N - 1)){
            if(matrix[i][M - 1] < 0){
                printf("%.15e]\n", matrix[i][M - 1]); 
            }else{
                printf("%.16e]\n", matrix[i][M - 1]); 
            }
        }else{
            if(matrix[i][M - 1] < 0){
                printf("%.15e;\n", matrix[i][M - 1]); 
            }else{
                printf("%.16e;\n", matrix[i][M - 1]); 
            }
        }
        
    }
    printf("\n"); 

    return; 
}

void print_matrix_i(int N, int M, int** matrix){
    printf("\n"); 
    for(int i = 0; i < N; i++){
        if(i == 0){
            printf("[");
        }else{
            printf(" ");
        }
        for(int j = 0; j < M - 1; j++){
            printf("%d ", matrix[i][j]); 
        }
        if(i == (N - 1)){
            printf("%d]\n", matrix[i][M - 1]); 
        }else{
            printf("%d;\n", matrix[i][M - 1]); 
        }
    }
    printf("\n"); 

    return; 
}

void read_measurement(int dim, const char* M_FILE, double* mu, double** S, double* T){
    FILE *m_file = fopen(M_FILE, "r");
    if (m_file == NULL) {
        fprintf(stderr, "Error: could not open file %s", M_FILE);
        exit(EXIT_FAILURE);
    }

    char line[256];
    char *token;
    int count = 0;

    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // mean vector
    token = strtok(line, " ");
    while (token != NULL && count < dim) { // read mean vector
        mu[count++] = strtod(token, NULL);
        token = strtok(NULL, " ");
    }
    count = 0;
    
    // Read covariance matrix
    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    for (int i = 0; i < dim; i++) { // read covariance matrix
        fgets(line, sizeof(line), m_file);
        token = strtok(line, " ");
        while (token != NULL && count < dim){
            S[i][count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }
        count = 0;
    }

    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // read T value
    *T = strtod(line, NULL);

    fclose(m_file);
    return;
}

void read_measurement_p(int dim, const char* M_FILE, double* mu, double** S, double* JD_0, double* T){
    FILE *m_file = fopen(M_FILE, "r");
    if (m_file == NULL) {
        fprintf(stderr, "Error: could not open file %s\n", M_FILE);
        exit(EXIT_FAILURE);
    }

    char line[256];
    char *token;
    char *saveptr; // thread-safe state for strtok_r
    int count;

    // Skip initial lines
    fgets(line, sizeof(line), m_file); // skip SPK ID line
    fgets(line, sizeof(line), m_file); // skip epoch_cal line
    fgets(line, sizeof(line), m_file); // epoch line

    // Extract JD_0 from line containing ':'
    char *p = strchr(line, ':');
    if (p == NULL) {
        fprintf(stderr, "Colon not found in line\n");
        exit(EXIT_FAILURE);
    }
    *JD_0 = strtod(p + 1, NULL);

    // Skip blank line and label line
    fgets(line, sizeof(line), m_file);
    fgets(line, sizeof(line), m_file);

    // Read mean vector
    fgets(line, sizeof(line), m_file);
    count = 0;
    token = strtok_r(line, " ", &saveptr);
    while (token != NULL && count < dim) {
        mu[count++] = strtod(token, NULL);
        token = strtok_r(NULL, " ", &saveptr);
    }

    // Skip blank line and label line
    fgets(line, sizeof(line), m_file);
    fgets(line, sizeof(line), m_file);

    // Read covariance matrix
    for (int i = 0; i < dim; i++) {
        fgets(line, sizeof(line), m_file);
        count = 0;
        char *row_saveptr = NULL;  // separate saveptr per row
        token = strtok_r(line, " ", &row_saveptr);
        while (token != NULL && count < dim) {
            S[i][count++] = strtod(token, NULL);
            token = strtok_r(NULL, " ", &row_saveptr);
        }
    }

    // Skip blank line, label line, then read T
    fgets(line, sizeof(line), m_file);
    fgets(line, sizeof(line), m_file);
    fgets(line, sizeof(line), m_file);
    *T = strtod(line, NULL);

    fclose(m_file);
}

void read_measurement_2_p(int dim, const char* M_FILE, double* mu, double* T){
    FILE *m_file = fopen(M_FILE, "r");
    if (m_file == NULL) {
        fprintf(stderr, "Error: could not open file %s\n", M_FILE);
        exit(EXIT_FAILURE);
    }

    char line[256];
    char *token;
    char *saveptr; // thread-safe state for strtok_r
    int count;

    // Skip initial lines
    fgets(line, sizeof(line), m_file); // skip label line

    // Read mean vector
    fgets(line, sizeof(line), m_file);
    count = 0;
    token = strtok_r(line, " ", &saveptr);
    while(token != NULL && count <= dim){
        if(count != dim){
            mu[count++] = strtod(token, NULL);
        }else{
            *T = strtod(token, NULL);
            count++;
        }
        token = strtok_r(NULL, " ", &saveptr);
    }

    fclose(m_file);
}

void read_measurement_g(int dim, int y_dim, const char* M_FILE, double* x_g, double* mu, double** S, double* T){
    FILE *m_file = fopen(M_FILE, "r");
    if (m_file == NULL) {
        fprintf(stderr, "Error: could not open file %s", M_FILE);
        exit(EXIT_FAILURE);
    }

    char line[256];
    char* token;
    int count = 0;

    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // mean vector
    token = strtok(line, " ");
    while (token != NULL && count < dim) { // read mean vector
        x_g[count++] = strtod(token, NULL);
        token = strtok(NULL, " ");
    }
    count = 0;

    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // mean vector
    token = strtok(line, " ");
    while (token != NULL && count < y_dim) { // read mean vector
        mu[count++] = strtod(token, NULL);
        token = strtok(NULL, " ");
    }
    count = 0;
    
    // Read covariance matrix
    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    for (int i = 0; i < y_dim; i++) { // read covariance matrix
        fgets(line, sizeof(line), m_file);
        token = strtok(line, " ");
        while (token != NULL && count < y_dim){
            S[i][count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }
        count = 0;
    }

    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // read T value
    *T = strtod(line, NULL);

    fclose(m_file);
    return;
}

void read_gbees(int N, int dim, const char* M_FILE, double** X, double* P){
    FILE *m_file = fopen(M_FILE, "r");
    if (m_file == NULL) {
        fprintf(stderr, "Error: could not open file %s", M_FILE);
        exit(EXIT_FAILURE);
    }

    char line[256];
    char *token;
    int count = 0;

    fgets(line, sizeof(line), m_file); // skip time line
    for(int i = 0; i < N; i++){
        fgets(line, sizeof(line), m_file); // read each row
        token = strtok(line, " ");
        P[i] = strtod(token, NULL); // read probability value
        token = strtok(NULL, " ");
        count = 0;
        while (token != NULL && count < dim){
            X[i][count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }
    }

    fclose(m_file);
    return;
}

char* JD_to_cal(double JD){

    double Z, F, A, B, C, D, E, alpha;
    JD += 0.5;
    Z = floor(JD);
    F = JD - Z;

    // Gregorian reform check
    if (Z < 2299161.0) {
        A = Z;
    } else {
        alpha = floor((Z - 1867216.25) / 36524.25);
        A = Z + 1 + alpha - floor(alpha / 4);
    }

    B = A + 1524;
    C = floor((B - 122.1) / 365.25);
    D = floor(365.25 * C);
    E = floor((B - D) / 30.6001);

    double day_frac = B - D - floor(30.6001 * E) + F;
    int day = (int)floor(day_frac);
    double frac_day = day_frac - day;

    int month = (E < 14) ? (int)(E - 1) : (int)(E - 13);
    int year = (month > 2) ? (int)(C - 4716) : (int)(C - 4715);

    // Convert fractional day into time (UTC-like)
    double total_seconds = frac_day * 86400.0;
    int hour = (int)(total_seconds / 3600.0);
    int minute = (int)((total_seconds - hour * 3600.0) / 60.0);
    int second = (int)(round(total_seconds - hour * 3600.0 - minute * 60.0));

    // Correct rounding overflow (e.g., 23:59:60)
    if (second == 60) {
        second = 0;
        minute++;
    }
    if (minute == 60) {
        minute = 0;
        hour++;
    }
    if (hour == 24) {
        hour = 0;
        day++;
        // Month/day overflow handling for completeness
        int days_in_month[] = { 0,31,28,31,30,31,30,31,31,30,31,30,31 };
        int leap = ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0));
        if (leap) days_in_month[2] = 29;
        if (day > days_in_month[month]) {
            day = 1;
            month++;
            if (month > 12) {
                month = 1;
                year++;
            }
        }
    }

    // Allocate formatted output string
    char *output = malloc(32 * sizeof(char));
    if (!output) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    snprintf(output, 32, "%04d-%02d-%02d %02d:%02d:%02d",
             year, month, day, hour, minute, second);

    return output;
}

double dot(double* vector_a, double* vector_b, int size){
    int i;
    double result = 0;
    for (i = 0; i < size; i++) {
        result += vector_a[i] * vector_b[i];
    }
    return result;
}

double** chol(int dim, double** S){
    double** L = init_m(dim, dim); 

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j <= i; j++){
            double sum = 0.0;
            for (int k = 0; k < j; k++){
                sum += L[i][k] * L[j][k];
            }
            if (i == j) {
                if (S[i][i] - sum < 0){
                    fprintf(stderr, "Matrix is not positive definite.\n");
                    return NULL; // Matrix not positive definite
                }
                L[i][j] = sqrt(S[i][i] - sum);
            } else {
                L[i][j] = (1.0 / L[j][j]) * (S[i][j] - sum);
            }
        }
    }
    return L;
}

int chol2(int dim, double* S, double* L){
    int flag = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', dim, S, dim);
    if(flag != 0){
        fprintf(stderr,
                "chol: Cholesky decomposition failed, flag = %d.\n",
                flag);
        free(S);
        return -2;
    }
    for (int i = 0; i < dim; ++i){
        for (int j = 0; j < dim; ++j){
            if (i >= j)
                L[i * dim + j] = S[i * dim + j];
            else
                L[i * dim + j] = 0.0;
        }
    }

    return flag; 
}

double** transpose(int N, int M, double** matrix){
    double** mT = init_m(M, N);

    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            mT[j][i] = matrix[i][j]; 
        }
    }

    return mT; 
}

void swap_rows(int dim, double** matrix, int row1, int row2){
    for (int col = 0; col < dim; ++col) {
        double temp = matrix[row1][col];
        matrix[row1][col] = matrix[row2][col];
        matrix[row2][col] = temp;
    }
}

void swap_vector_rows(int dim, double** vector, int row1, int row2){
    for (int col = 0; col < dim; ++col) {
        double temp = vector[row1][col];
        vector[row1][col] = vector[row2][col];
        vector[row2][col] = temp;
    }
}

double** inv(int dim, double** mat){
    double** inverse = init_m(dim, dim);
    double** matrix = init_m(dim, dim);

    // Initialize the inverse matrix as the identity matrix
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            inverse[i][j] = (i == j) ? 1.0 : 0.0;
            matrix[i][j] = mat[i][j];
        }
    }

    // Perform Gaussian elimination
    for (int i = 0; i < dim; ++i) {
        // Find the row with the maximum element in the current column
        double max_val = fabs(matrix[i][i]);
        int max_row = i;
        for (int k = i + 1; k < dim; ++k) {
            if (fabs(matrix[k][i]) > max_val) {
                max_val = fabs(matrix[k][i]);
                max_row = k;
            }
        }

        // If the maximum value is zero, the matrix is singular
        if (fabs(max_val) < DBL_MIN){
            // print_matrix(dim, dim, mat); 
            printf("Matrix is singular and cannot be inverted.\n");
            return NULL;
        }

        // Swap the current row with the row of the maximum element
        if (max_row != i) {
            swap_rows(dim, matrix, i, max_row);
            swap_vector_rows(dim, inverse, i, max_row);
        }

        // Normalize the current row
        double diagonal_val = matrix[i][i];
        for (int col = 0; col < dim; ++col) {
            matrix[i][col] /= diagonal_val;
            inverse[i][col] /= diagonal_val;
        }

        // Eliminate the current column in all other rows
        for (int k = 0; k < dim; ++k) {
            if (k != i) {
                double factor = matrix[k][i];
                for (int col = 0; col < dim; ++col) {
                    matrix[k][col] -= factor * matrix[i][col];
                    inverse[k][col] -= factor * inverse[i][col];
                }
            }
        }
    }
    return inverse; 
}

double generate_standard_normal(void){
    double u1 = rand() / (RAND_MAX + 1.0);
    double u2 = rand() / (RAND_MAX + 1.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}


double generate_standard_normal_p(unsigned int* seed){
    double u1 = rand_r(seed) / (RAND_MAX + 1.0);
    double u2 = rand_r(seed) / (RAND_MAX + 1.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

double* mvnrnd(int N, int M, const double* mu, double** S){
    double** L = chol(M, S); 
    double* samples = init_v(N * M); 
    for (int n = 0; n < N; n++) {
        double z[M];
        // Generate standard normal vector
        for (int m = 0; m < M; m++) {
            z[m] = generate_standard_normal();
        }
        // Transform the vector
        for (int m = 0; m < M; m++) {
            samples[n * M + m] = mu[m];
            for (int j = 0; j <= m; j++) {
                samples[n * M + m] += L[m][j] * z[j];
            }
        }
    }

    return samples;
}

double* mvnrnd_p(int N, int M, const double* mu, double** S, unsigned int* seed){
    double** L = chol(M, S); 
    if(!L){
        #pragma omp critical
        {
            printf("Cholesky failed for this covariance matrix!\n");
        }
        return NULL;
    }
    double* samples = init_v(N * M); 
    for (int n = 0; n < N; n++) {
        double z[M];
        // Generate standard normal vector
        for (int m = 0; m < M; m++) {
            z[m] = generate_standard_normal_p(seed);
        }
        // Transform the vector
        for (int m = 0; m < M; m++) {
            samples[n * M + m] = mu[m];
            for (int j = 0; j <= m; j++) {
                samples[n * M + m] += L[m][j] * z[j];
            }
        }
    }
    for(int i=0; i<M; i++) free(L[i]);
    free(L);

    return samples;
}

void RK4(void (*f)(double t, double* x, double* xp, double* coef), int dim, double* x, double* coef, double t, double dt){

  double* f1 = init_v(dim); 
  double* f2 = init_v(dim); 
  double* f3 = init_v(dim); 
  double* f4 = init_v(dim); 
  double* x0 = init_v(dim); 
  for(int i = 0; i < dim; i++) x0[i] = x[i]; 

  f(t, x, f1, coef);            for(int i = 0; i < dim; i++) x[i] = x0[i] + (dt / 2) * f1[i]; 
  f(t + (dt / 2), x, f2, coef); for(int i = 0; i < dim; i++) x[i] = x0[i] + (dt / 2) * f2[i]; 
  f(t + (dt / 2), x, f3, coef); for(int i = 0; i < dim; i++) x[i] = x0[i] + dt * f3[i]; 
  f(t + dt, x, f4, coef);       for(int i = 0; i < dim; i++) x[i] = x0[i] + dt * ((1.0 / 6.0) * f1[i] + (1.0 / 3.0) * f2[i] + (1.0 / 3.0) * f3[i] + (1.0 / 6.0) * f4[i]);

  free(f1); 
  free(f2); 
  free(f3); 
  free(f4); 
  free(x0);

  return; 
}

double gauss_probability(int dim, double* x, double* mu, double** S){ 
    double diff[dim];
    for(int i=0; i<dim; i++) diff[i] = x[i] - mu[i];

    double** M_inv = inv(dim, S);
    double* M_inv_x = mvtimes(dim, dim, M_inv, diff);
    double dot_prod = dot(diff, M_inv_x, dim);

    for(int i=0; i<dim; i++) free(M_inv[i]);
    free(M_inv);
    free(M_inv_x);

    return exp(-0.5*dot_prod);
}

double* mean(int N, int M, double* X, double* p){
    double* m = init_v(M);

    if(p == NULL){
        p = init_v(N);
        for(int i = 0; i < N; i++){
            p[i] = 1.0 / N; 
        }
    }

    for(int i = 0; i < M; i++){
        m[i] = 0; 
        for(int j = 0; j < N; j++){
            m[i] += p[j] * X[j * M + i]; 
        }
    }

    return m; 
}

double** cov(int N, int M, double* X, double* p){
    double** S = init_m(M, M);

    if (p == NULL) {
        p = init_v(N);
        for (int i = 0; i < N; i++)
            p[i] = 1.0 / N;
    }

    /* ---- NaN check on X ---- */
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < M; i++) {
            if (isnan(X[k * M + i])) {
                fprintf(stderr,
                        "ERROR: NaN detected in X[%d][%d]\n", k, i);
                abort();
            }
        }
    }

    /* ---- NaN check on p ---- */
    for (int k = 0; k < N; k++) {
        if (isnan(p[k])) {
            fprintf(stderr,
                    "ERROR: NaN detected in p[%d]\n", k);
            abort();
        }
    }

    double* m = mean(N, M, X, p);

    /* ---- NaN check on mean ---- */
    for (int i = 0; i < M; i++) {
        if (isnan(m[i])) {
            fprintf(stderr,
                    "ERROR: NaN detected in mean[%d]\n", i);
            abort();
        }
    }

    /* ---- Covariance ---- */
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            S[i][j] = 0.0;
            for (int k = 0; k < N; k++) {
                double di = X[k * M + i] - m[i];
                double dj = X[k * M + j] - m[j];

                if (isnan(di) || isnan(dj)) {
                    fprintf(stderr,
                        "ERROR: NaN during covariance at "
                        "X[%d][%d], X[%d][%d]\n",
                        k, i, k, j);
                    abort();
                }

                S[i][j] += p[k] * di * dj;
            }

            if (isnan(S[i][j])) {
                fprintf(stderr,
                        "ERROR: NaN in covariance S[%d][%d]\n", i, j);
                abort();
            }
        }
    }

    /* ---- Enforce symmetry ---- */
    const double SYM_TOL = 1e-12;

    for (int i = 0; i < M; i++) {
        for (int j = i + 1; j < M; j++) {
            double a = S[i][j];
            double b = S[j][i];

            if (fabs(a - b) > SYM_TOL) {
                fprintf(stderr,
                        "WARNING: Covariance asymmetry S[%d][%d]=%.3e, "
                        "S[%d][%d]=%.3e (symmetrizing)\n",
                        i, j, a, j, i, b);
            }

            double s = 0.5 * (a + b);
            S[i][j] = s;
            S[j][i] = s;
        }
    }

    free(m);
    return S;
}

double mahalanobis_distance(double *xi, double *xj, double **S_inv, int dim){
    double* diff = init_v(dim); 
    for (int k = 0; k < dim; ++k)
        diff[k] = xi[k] - xj[k];

    // temp = S_inv * diff
    double *temp = (double *)malloc(dim * sizeof(double));
    for (int i = 0; i < dim; ++i) {
        temp[i] = 0.0;
        for (int j = 0; j < dim; ++j)
            temp[i] += S_inv[i][j] * diff[j];
    }

    // result = diff^T * temp
    double result = 0.0;
    for (int i = 0; i < dim; ++i)
        result += diff[i] * temp[i];

    free(diff);
    free(temp);
    return result;
}

double* cholesky_upper(int dim, double** S){
    double* R = (double*)calloc(dim*dim, sizeof(double));
    if(!R) return NULL;
    for(int i=0;i<dim;++i){
        for(int j=i;j<dim;++j){
            double sum = S[i][j];
            for(int k=0;k<i;++k) sum -= R[k*dim + i] * R[k*dim + j];
            if(i==j){
                if(sum <= 0.0){ free(R); return NULL; } /* not PD */
                R[i*dim + i] = sqrt(sum);
            } else {
                R[i*dim + j] = sum / R[i*dim + i];
            }
        }
    }
    return R;
}

void solve_Rt(int dim, const double* R, const double* b, double* out){
    /* Solve (R^T) x = b  -> lower triangular solve */
    for(int i=0;i<dim;++i){
        double s = b[i];
        for(int k=0;k<i;++k) s -= R[k*dim + i] * out[k]; /* R[k,i] stored at R[k*dim + i] */
        out[i] = s / R[i*dim + i];
    }
}

/* compute squared Euclidean distance between two length-d vectors */
double sq_euclid(int dim, const double* a, const double* b){
    double s=0.0;
    for(int k=0;k<dim;++k){
        double t = a[k] - b[k];
        s += t*t;
    }
    return s;
}

double* flatten_X(int N, int dim, double** X){
    double* out = (double*)malloc(sizeof(double)*N*dim);
    if(!out) return NULL;
    for(int i=0;i<N;++i){
        memcpy(out + i*dim, X[i], sizeof(double)*dim);
    }
    return out;
}

double calc_HZ_fast(int N, int dim, double* X, double* p, double* mu, double** S){
    if(N <= 0 || dim <= 0) return 0.0;

    /* compute beta and frequently-used constants once */
    double beta = (1.0 / sqrt(2.0)) * pow(((N * (2.0 * dim + 1.0)) / 4.0), 1.0 / (dim + 4.0));
    double beta2 = beta * beta;
    double c1 = beta2 * 0.5; /* for HZ1: exp(-c1 * d) where d is Mahalanobis quadratic */
    double c2 = beta2 / (2.0 * (1.0 + beta2)); /* for HZ2 exponent */
    double two_times = 2.0 * pow((1.0 + beta2), -((double)dim)/2.0);
    double last_term = N * pow((1.0 + 2.0 * beta2), -((double)dim)/2.0);

    /* compute mean and covariance (user-supplied) */
    double* mu_t = mean(N, dim, X, p); /* length dim, must be freed by us */
    double** S_t = cov(N, dim, X, p);  /* dim x dim double** (row pointers), must be freed by us */
    for(int i = 0; i < dim; i++){
        mu[i] = mu_t[i]; 
        for(int j = 0; j < dim; j++){
            S[i][j] = S_t[i][j]; 
        }
    }
    free(mu_t);
    for(int i = 0; i < dim; i++){
        free(S_t[i]); 
    }
    free(S_t); 

    /* Cholesky of S -> R (upper triangular) */
    double* R = cholesky_upper(dim, S);
    if(!R){
        /* fallback: compute S inverse and do original way (not implemented here) */
        /* free and return error-like value */
        /* free S row pointers if needed according to your cov() alloc scheme */
        return NAN;
    }

    /* Transform points: z_i = R^{-T} x_i  (i.e. solve R^T z = x) */
    double* Z = (double*)malloc(sizeof(double)*N*dim);
    if(!Z){ free(X); free(R); return NAN; }

    double* temp = (double*)malloc(sizeof(double)*dim);
    for(int i=0;i<N;++i){
        solve_Rt(dim, R, X + i*dim, temp);
        memcpy(Z + i*dim, temp, sizeof(double)*dim);
    }
    /* transform mu as well */
    double* zmu = (double*)malloc(sizeof(double)*dim);
    solve_Rt(dim, R, mu, zmu);

    free(temp);
    
    /* Compute HZ1 using symmetric pair sum to cut roughly half computations.
       Sum_{i,j} exp(-c1 * ||z_i - z_j||^2)  = N (i==j terms) + 2 * sum_{i<j} exp(...)
    */
    double HZ1 = 0.0;
    /* diagonal terms contribute exp(0) = 1 each -> initial N */
    HZ1 += (double)N;

    /* accumulate off-diagonals */
    double off = 0.0;
    /* Optional: enable OpenMP parallelization here for large N by compiling with -fopenmp
       and including #include <omp.h>. Then use pragma below.
     */
    #pragma omp parallel for reduction(+:off) schedule(static)
    for(int i=0;i<N;++i){
        const double* zi = Z + i*dim;
        for(int j=i+1;j<N;++j){
            const double* zj = Z + j*dim;
            double d2 = sq_euclid(dim, zi, zj); /* squared Euclidean in transformed space == Mahalanobis */
            off += exp(-c1 * d2);
        }
    }
    HZ1 += 2.0 * off;
    HZ1 /= (double)N;

    /* Compute HZ2: sum_i exp(-c2 * ||z_i - zmu||^2) then multiply by two_times */
    double HZ2 = 0.0;
    #pragma omp parallel for reduction(+:HZ2) schedule(static)
    for(int i=0;i<N;++i){
        const double* zi = Z + i*dim;
        double d2 = sq_euclid(dim, zi, zmu);
        HZ2 += exp(-c2 * d2);
    }
    HZ2 *= two_times;

    /* final HZ */
    double HZ = HZ1 - HZ2 + last_term;

    /* cleanup */
    free(Z);
    free(zmu);
    free(R);

    return HZ;
}

int logspace_fill(double a, double b, int n, double* out){
    if (n <= 0 || out == NULL) return -1;
    if (n == 1) {
        out[0] = pow(10.0, a);
        return 0;
    }

    double da = (b - a) / (double)(n - 1);
    for (int i = 0; i < n; ++i) {
        double expv = a + da * i;                /* exponent in base-10 */
        out[i] = pow(10.0, expv);               /* 10^(a + i*da) */
    }
    return 0;
}

double* logspace(double a, double b, int n){
    if (n <= 0) return NULL;
    double* arr = init_v(n); 
    if(logspace_fill(a, b, n, arr) != 0) {free(arr); return NULL;}
    return arr;
}

int count_txt_files(const char *path){
    DIR *dir;
    struct dirent *entry;
    int count = 0;

    dir = opendir(path);
    if (!dir) {
        perror("opendir");
        return -1;
    }

    while ((entry = readdir(dir)) != NULL) {
        // Skip directories and hidden files
        if (entry->d_type == DT_DIR)
            continue;

        // Find the file extension
        const char *ext = strrchr(entry->d_name, '.');
        if (ext && strcmp(ext, ".txt") == 0)
            count++;
    }

    closedir(dir);
    return count;
}

void measurement_update_gbees(void (*h)(double t, double* x, double* y, double* coef), double** X, double* P, double* y, double* coef, double** Sy, double t, double thresh, int N, int dim, int y_dim, double*** Xp, double** Pp, int* Np){

    // likelihood
    double* l = init_v(N);
    double* y_i = init_v(y_dim); 
    double sum = 0.0; 
    for(int i = 0; i < N; i++){
        h(t, X[i], y_i, coef); 
        l[i] = gauss_probability(y_dim, y_i, y, Sy);
        P[i] *= l[i];
        sum += P[i];
    }

    // renormalize weights
    *Np = 0;
    for(int i = 0; i < N; i++){
        P[i] /= sum; 
        if(P[i] > thresh) (*Np)++;
    }

    *Xp = init_m(*Np, dim);
    *Pp = init_v(*Np);
    sum = 0.0; 
    int count = 0; 
    for(int i = 0; i < N; i++){
        if(P[i] > thresh){
            for(int j = 0; j < dim; j++){
                (*Xp)[count][j] = X[i][j]; 
            }
            (*Pp)[count] = P[i];
            sum += (*Pp)[count++];
        }
    }
    for(int i = 0; i < *Np; i++){
        (*Pp)[i] /= sum; 
    }

    free(l); 
    free(y_i); 

    return;
}

double wrap_angle(double a){
    return atan2(sin(a), cos(a));
}

double measurement_residual(double a, double b, int index){
    double d = a - b;

    // y[1] is theta
    if (index == 1) {
        return wrap_angle(d);
    }

    return d;
}

double norm_inf(double* g, int p){
    double norm = 0.0; 
    for(int i = 0; i < p; i++){
        double val = fabs(g[i]);
        if(val > norm){
            norm = val;
        }
    }
    return norm; 
}

double norm_frob(double* g, int p){
    double norm = 0.0; 
    for(int i = 0; i < p; i++){
        norm += (fabs(g[i]) * fabs(g[i]));
    }
    return pow(norm, 0.5); 
}

int expm(double* A, double* E, int n){

    const double theta13 = 5.371920351148152;
    const double b[] = {
        64764752532480000.0,
        32382376266240000.0,
        7771770303897600.0,
        1187353796428800.0,
        129060195264000.0,
        10559470521600.0,
        670442572800.0,
        33522128640.0,
        1323241920.0,
        40840800.0,
        960960.0,
        16380.0,
        182.0,
        1.0
    };

    int nn = n * n;

    double* As = init_v(nn);
    double* A2 = init_v(nn); 
    double* A4 = init_v(nn);
    double* A6 = init_v(nn);
    double* U = init_v(nn);
    double* V = init_v(nn);
    double* M1 = init_v(nn);
    double* M2 = init_v(nn);

    if (!As || !A2 || !A4 || !A6 || !U || !V || !M1 || !M2)
        return -1;

    memcpy(As, A, nn * sizeof(double));

    /* 1-norm */
    double norm = 0.0;

    for(int j = 0; j < n; j++){
        double s = 0.0;
        for(int i = 0; i < n; i++) s += fabs(As[i*n+j]);
        if(s > norm) norm = s;
    }

    int s = 0;
    if(norm > theta13) s = (int)ceil(log2(norm / theta13));
    double scale = ldexp(1.0, -s);
    for(int i = 0; i < nn; i++) As[i] *= scale;

    /* A2 = A^2 */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, As, n, As, n, 0.0, A2, n);

    /* A4 = A2^2 */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A2, n, A2, n, 0.0, A4, n);

    /* A6 = A4*A2 */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A4, n, A2, n, 0.0, A6, n);

    /*
     * U = A * [A6*(b13*A6+b11*A4+b9*A2)
     *          + b7*A6+b5*A4+b3*A2+b1*I]
     */
    for (int i = 0; i < nn; i++) M1[i] = b[13]*A6[i] + b[11]*A4[i] + b[9]*A2[i];

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A6, n, M1, n, 0.0, M2, n);

    for (int i = 0; i < nn; i++) M2[i] += b[7]*A6[i] + b[5]*A4[i] + b[3]*A2[i];

    for (int i = 0; i < n; i++) M2[i*n+i] += b[1];

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, As, n, M2, n, 0.0, U, n);

    /*
     * V = A6*(b12*A6+b10*A4+b8*A2)
     *     + b6*A6+b4*A4+b2*A2+b0*I
     */
    for(int i = 0; i < nn; i++) M1[i] = b[12]*A6[i] + b[10]*A4[i] + b[8]*A2[i];

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A6, n, M1, n, 0.0, V, n);

    for(int i = 0; i < nn; i++) V[i] += b[6]*A6[i] + b[4]*A4[i] + b[2]*A2[i];

    for(int i = 0; i < n; i++) V[i*n+i] += b[0];

    /*
     * Solve:
     *
     *     (V-U) E = (V+U)
     */

    for(int i = 0; i < nn; i++){
        M1[i] = V[i] - U[i];
        M2[i] = V[i] + U[i];
    }

    int *ipiv = malloc(n*sizeof(int));

    if(!ipiv){
        free(As); free(A2); free(A4); free(A6);
        free(U); free(V); free(M1); free(M2);
        return -1;
    }

    int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, n, M1, n, ipiv, M2, n);
    free(ipiv);

    if(info != 0){
        free(As); free(A2); free(A4); free(A6);
        free(U); free(V); free(M1); free(M2);
        return info;
    }

    memcpy(E, M2, nn * sizeof(double));

    /* Squaring */
    for(int k = 0; k < s; k++){

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, E, n, E, n, 0.0, M1, n);

        memcpy(E, M1, nn*sizeof(double));
    }

    free(As);
    free(A2);
    free(A4);
    free(A6);
    free(U);
    free(V);
    free(M1);
    free(M2);

    return 0;
}

double** dt_process_noise(void (*J)(double t, double* x, double** F), double t, double** Q, double* x, int dim){

    double** M1 = init_m(dim, dim);
    double** M3 = init_m(dim, dim);
    double** F  = init_m(dim, dim);

    J(t, x, F);

    /* Van Loan blocks:
     * A = [ -F      G Q G^T
     *        0       F^T  ]
     */
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            M1[i][j] = -F[i][j];
            M3[i][j] =  F[j][i];
        }
    }

    int dof = dim / 2;
    double** G  = init_m(dim, dof);
    double** Gt = init_m(dof, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dof; j++) {
            if (i < dof) {
                G[i][j] = 0.0;
                Gt[j][i] = 0.0;
            }else if((i - dof) == j){
                G[i][j] = 1.0;
                Gt[j][i] = 1.0;
            }else{
                G[i][j] = 0.0;
                Gt[j][i] = 0.0;
            }
        }
    }

    /* G Q G^T */
    double** M2_a = mmtimes(dof, dof, dim, Q, Gt);
    double** M2   = mmtimes(dim, dof, dim, G, M2_a);

    int n2 = 2 * dim;
    double* A    = init_v(n2 * n2);
    double* expA = init_v(n2 * n2);

    /* Construct Van Loan matrix */
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {

            A[i * n2 + j] =
                M1[i][j];

            A[i * n2 + (j + dim)] =
                M2[i][j];

            A[(i + dim) * n2 + j] =
                0.0;

            A[(i + dim) * n2 + (j + dim)] =
                M3[i][j];
        }
    }

    /* exp(A) */
    int flag = expm(A, expA, n2);
    if(flag != 0){
        fprintf(stderr, "expm failed: %d\n", flag);
        return NULL;
    }

    /* Extract blocks from exp(A) */
    double** U   = init_m(dim, dim);
    double** PhiT = init_m(dim, dim);

    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            U[i][j] = expA[i * n2 + (j + dim)];
            PhiT[i][j] = expA[(i + dim) * n2 + (j + dim)];
        }
    }

    /*
     * Q_k = Phi * U
     *
     * Since PhiT = Phi^T,
     * Phi = PhiT^T.
     */
    double** Phi = transpose(dim, dim, PhiT);
    double** Qk  = mmtimes(dim, dim, dim, Phi, U);

    return Qk;
}

int parse_gbees(const char* filename, int dim, double** P, double*** X, int* N){

    FILE* fid = fopen(filename, "rb");
    if (fid == NULL) {
        fprintf(stderr, "Could not open file: %s\n", filename);
        return -1;
    }

    double t;
    if (fread(&t, sizeof(double), 1, fid) != 1) {
        fprintf(stderr, "Error reading t.\n");
        fclose(fid);
        return -1;
    }

    uint32_t n_uint;
    if (fread(&n_uint, sizeof(uint32_t), 1, fid) != 1) {
        fprintf(stderr, "Error reading n.\n");
        fclose(fid);
        return -1;
    }

    *N = (int)n_uint;
    *P = init_v(*N);
    *X = init_m(*N, dim);
    for (int i = 0; i < *N; i++) {

        float value;

        /* P */
        if (fread(&value, sizeof(float), 1, fid) != 1) {
            fprintf(stderr, "Error reading P[%d].\n", i);
            goto cleanup;
        }

        (*P)[i] = (double)value;

        /* x */
        for (int j = 0; j < dim; j++) {

            if (fread(&value, sizeof(float), 1, fid) != 1) {
                fprintf(stderr, "Error reading x[%d][%d].\n", i, j);
                goto cleanup;
            }

            (*X)[i][j] = (double)value;
        }
    }

    fclose(fid);
    return 0;


cleanup:

    for (int i = 0; i < *N; i++) {
        free((*X)[i]);
    }

    free(*X);
    free(*P);

    fclose(fid);

    return -1;
}

void weighted_mean_cov(double** X, double* P, int N, int dim, double* mu, double** S){

    double Psum = 0.0;
    for(int i = 0; i < N; i++){
        Psum += P[i];
    }

    /* Initialize mean */
    for(int j = 0; j < dim; j++){
        mu[j] = 0.0;
    }

    /* Weighted mean */
    for(int i = 0; i < N; i++){
        for(int j = 0; j < dim; j++){
            mu[j] += P[i] * X[i][j];
        }
    }

    /* Initialize covariance */
    for(int j = 0; j < dim; j++){
        for(int k = 0; k < dim; k++){
            S[j][k] = 0.0;
        }
    }

    /* Weighted covariance */
    for(int i = 0; i < N; i++){
        double w = P[i] / Psum;
        for(int j = 0; j < dim; j++){
            double dx_j = X[i][j] - mu[j];
            for(int k = 0; k < dim; k++){
                double dx_k = X[i][k] - mu[k];
                S[j][k] += w * dx_j * dx_k;
            }
        }
    }

    /* Forcing symmetry */
    for(int i = 0; i < dim; i++){
        for(int j = 0; j <= i; j++){
            double sym_val = 0.5 * (S[i][j] + S[j][i]);
            S[i][j] = sym_val;
            S[j][i] = sym_val;
        }
    }

    return; 
}

void weighted_mean_cov_2(double* X, double* P, int N, int dim, double* mu, double** S){

    double Psum = 0.0;
    for(int i = 0; i < N; i++){
        Psum += P[i];
    }

    /* Initialize mean */
    for(int j = 0; j < dim; j++){
        mu[j] = 0.0;
    }

    /* Weighted mean */
    for(int i = 0; i < N; i++){
        for(int j = 0; j < dim; j++){
            mu[j] += P[i] * X[i * dim + j];
        }
    }

    /* Initialize covariance */
    for(int j = 0; j < dim; j++){
        for(int k = 0; k < dim; k++){
            S[j][k] = 0.0;
        }
    }

    /* Weighted covariance */
    for(int i = 0; i < N; i++){
        double w = P[i] / Psum;
        for(int j = 0; j < dim; j++){
            double dx_j = X[i * dim + j] - mu[j];
            for(int k = 0; k < dim; k++){
                double dx_k = X[i * dim + k] - mu[k];
                S[j][k] += w * dx_j * dx_k;
            }
        }
    }

    /* Forcing symmetry */
    for(int i = 0; i < dim; i++){
        for(int j = 0; j <= i; j++){
            double sym_val = 0.5 * (S[i][j] + S[j][i]);
            S[i][j] = sym_val;
            S[j][i] = sym_val;
        }
    }

    return; 
}

void write_pdf(char* FILE_NAME, double t, int N, int dim, double* X, double* P){

    FILE* f = fopen(FILE_NAME, "wb");
    fwrite(&t, sizeof(double), 1, f);
    uint32_t n_uint = (uint32_t)N;
    fwrite(&n_uint, sizeof(uint32_t), 1, f);
    for (int i = 0; i < N; i++) {
        fwrite(&P[i], sizeof(double), 1, f);
        for (int j = 0; j < dim; j++) {
            fwrite(&X[i * dim + j], sizeof(double), 1, f);
        }
    }
    fclose(f);
    return;
}

// void write_pdf(char* file_str, int N, int dim, double* X, double* P){

//     char FILE_NAME[256];  
//     snprintf(FILE_NAME, sizeof(FILE_NAME), "%s", file_str);
//     FILE* f = fopen(FILE_NAME, "w");
//     fprintf(f, "P "); 
//     for(int i = 1; i < dim; i++){
//         fprintf(f, "x%d ", i); 
//     }
//     fprintf(f, "x%d\n", dim); 
//     for(int i = 0; i < N; i++){
//         fprintf(f, "%.15e ", P[i]);
//         for(int j = 0; j < dim - 1; j++){
//             fprintf(f, "%.15e ", X[i * dim + j]); 
//         }
//         fprintf(f, "%.15e\n", X[i * dim + dim - 1]); 
//     }

//     return; 
// }

void write_pdf_2(char* file_str, int N, int dim, double* X, double* P, double* mu, double** S){

    char FILE_NAME[256];  
    snprintf(FILE_NAME, sizeof(FILE_NAME), "%s", file_str);
    FILE* f = fopen(FILE_NAME, "w");
    fprintf(f, "mu\n"); 
    for(int i = 0; i < dim; i++){
        fprintf(f, "%.15f ", mu[i]); 
    }
    fprintf(f, "\nS\n"); 
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            fprintf(f, "%.15f ", S[i][j]);
        }
        fprintf(f, "\n"); 
    }
    fprintf(f, "\nP "); 
    for(int i = 1; i < dim; i++){
        fprintf(f, "x%d ", i); 
    }
    fprintf(f, "x%d\n", dim); 
    for(int i = 0; i < N; i++){
        fprintf(f, "%.15e ", P[i]);
        for(int j = 0; j < dim - 1; j++){
            fprintf(f, "%.15e ", X[i * dim + j]); 
        }
        fprintf(f, "%.15e\n", X[i * dim + dim - 1]); 
    }

    return; 
}

double det(int dim, const double* M){
    double det = 1.0;
    for(int i = 0; i < dim; ++i){
        det *= M[i * dim + i]; // Multiply diagonal elements
    }
    return det;
}

void calc_b1_b2(int N, int dim, double* X, double* mu, double** S, double* b1, double* b2){

    double* mu_t = mean(N, dim, X, NULL);
    double** S_p = cov(N, dim, X, NULL);
    double* S_t = init_v(dim * dim); 
    for(int i = 0; i < dim; i++){
        mu[i] = mu_t[i]; 
        for(int j = 0; j < dim; j++){
            S[i][j] = S_p[i][j]; 
            S_t[i * dim + j] = S_p[i][j]; 
        }
    }

    double* Xc = init_v(N * dim);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            Xc[i * dim + j] = X[i * dim + j] - mu_t[j];
        }
    }


    double* A = init_v(dim * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            A[j * N + i] = Xc[i * dim + j];
        }
    }

    /*
     * -------------------------------------------------------
     *
     * dposv overwrites A and B.
     *
     * S * A = Xc'
     * 
     * -------------------------------------------------------
     */
    int info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', dim, N, S_t, dim, A, N);
    if(info != 0){
        fprintf(stderr, "LAPACKE_dposv failed, info = %d\n", info);
        free(mu_t);
        free(S_t);
        free(Xc);
        free(A);
        *b1 = NAN;
        *b2 = NAN;
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
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, dim, 1.0, Xc, dim, A, N, 0.0, Q, N);

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

    *b1 = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            const double q = Q[i*N + j];
            *b1 += q * q * q;
        }
    }
    *b1 /= (N * N); 

    /*
     * -------------------------------------------------------
     * MATLAB:
     *
     * bt2 = sum(W .* (q_diag.^2));
     * -------------------------------------------------------
     */

    *b2 = 0.0;
    for (int i = 0; i < N; i++){
        const double qii = Q[i * N + i];
        *b2 += qii * qii;
    }
    *b2 /= N; 

    free(mu_t);
    free(S_t);
    free(Xc);
    free(A);
    free(Q);

    return;
}

double calc_b1(int N, int dim, double* X, double* mu, double** S){

    double* mu_t = mean(N, dim, X, NULL);
    double** S_p = cov(N, dim, X, NULL);
    double* S_t = init_v(dim * dim); 
    for(int i = 0; i < dim; i++){
        mu[i] = mu_t[i]; 
        for(int j = 0; j < dim; j++){
            S[i][j] = S_p[i][j]; 
            S_t[i * dim + j] = S_p[i][j]; 
        }
    }

    double* Xc = init_v(N * dim);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            Xc[i * dim + j] = X[i * dim + j] - mu_t[j];
        }
    }


    double* A = init_v(dim * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            A[j * N + i] = Xc[i * dim + j];
        }
    }

    /*
     * -------------------------------------------------------
     *
     * dposv overwrites A and B.
     *
     * S * A = Xc'
     * 
     * -------------------------------------------------------
     */
    int info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', dim, N, S_t, dim, A, N);
    if(info != 0){
        fprintf(stderr, "LAPACKE_dposv failed, info = %d\n", info);
        free(mu_t);
        free(S_t);
        free(Xc);
        free(A);
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
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, dim, 1.0, Xc, dim, A, N, 0.0, Q, N);

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

    double b1 = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            const double q = Q[i*N + j];
            b1 += q * q * q;
        }
    }
    b1 /= (N * N); 

    free(mu_t);
    free(S_t);
    free(Xc);
    free(A);
    free(Q);

    return b1;
}


double calc_b2(int N, int dim, double* X, double* mu, double** S){

    double* mu_t = mean(N, dim, X, NULL);
    double** S_p = cov(N, dim, X, NULL);
    double* S_t = init_v(dim * dim); 
    for(int i = 0; i < dim; i++){
        mu[i] = mu_t[i]; 
        for(int j = 0; j < dim; j++){
            S[i][j] = S_p[i][j]; 
            S_t[i * dim + j] = S_p[i][j]; 
        }
    }

    double* Xc = init_v(N * dim);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            Xc[i * dim + j] = X[i * dim + j] - mu_t[j];
        }
    }


    double* A = init_v(dim * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            A[j * N + i] = Xc[i * dim + j];
        }
    }

    /*
     * -------------------------------------------------------
     *
     * dposv overwrites A and B.
     *
     * S * A = Xc'
     * 
     * -------------------------------------------------------
     */
    int info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', dim, N, S_t, dim, A, N);
    if(info != 0){
        fprintf(stderr, "LAPACKE_dposv failed, info = %d\n", info);
        free(mu_t);
        free(S_t);
        free(Xc);
        free(A);
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
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, dim, 1.0, Xc, dim, A, N, 0.0, Q, N);

    /*
     * -------------------------------------------------------
     *
     * q_diag = diag(Q);
     * -------------------------------------------------------
     */

    /*
     * -------------------------------------------------------
     * MATLAB:
     *
     * bt2 = sum(W .* (q_diag.^2));
     * -------------------------------------------------------
     */

    double b2 = 0.0;
    for (int i = 0; i < N; i++){
        const double qii = Q[i * N + i];
        b2 += qii * qii;
    }
    b2 /= N; 

    free(mu_t);
    free(S_t);
    free(Xc);
    free(A);
    free(Q);

    return b2;
}

bool check_b1_b2(double b1, double b2, int N, int dim, double z1_s, double z2_s){

    double z1 = ((dim + 1.0) * (N + 1.0) * (N + 3.0) / (6.0 * ((N + 1.0) * (dim + 1.0) - 6.0))) * b1;
    double z2 = (b2 - dim * (dim + 2.0)) / sqrt(8.0 * dim * (dim + 2.0) / N);

    return (z1 <= z1_s) && (z2 <= z2_s) && (z2 >= -z2_s);
}

bool check_b1(double b1, int N, int dim, double z1_s){

    double z1 = ((dim + 1.0) * (N + 1.0) * (N + 3.0) / (6.0 * ((N + 1.0) * (dim + 1.0) - 6.0))) * b1;

    return (z1 <= z1_s);
}

bool check_b2(double b2, int N, int dim, double z2_s){

    double z2 = (b2 - dim * (dim + 2.0)) / sqrt(8.0 * dim * (dim + 2.0) / N);

    return (z2 <= z2_s) && (z2 >= -z2_s);
}

double frob_norm(int n, int m, double** A){
    double sum = 0.0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            sum += A[i][j] * A[i][j];
        }
    }
    return sqrt(sum);
}