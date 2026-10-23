#include "../../src/utils.h"
#include "../../src/models.h"
#include "../../src/cukf.h"
#include "../../src/pme.h"

int main(int argc, char* argv[]){
    
    // inputs
    fn f = PCR3BP;
    int dim = 4;     
    fn h = R;
    const char* M_MODEL = "PCR3BP";
    int nm = 3; 
    int CUT_M = 6; 

    // parameters
    int dim     = 4;              
    double* mu  = init_v(dim); 
    double** S  = init_m(dim, dim);
    double T    = 0; 
    
    // directories
    char M_DIR[256];  
    snprintf(M_DIR, sizeof(M_DIR), "../../measurements/%s_hybrid", M_MODEL);
    char R_DIR[256];  
    snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s", M_MODEL);
    mkdir(R_DIR, 0777);
    snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s/cukf_gbees", M_MODEL);
    mkdir(R_DIR, 0777);
    char M_FILE[256];  
    snprintf(M_FILE, sizeof(M_FILE), "%s/measurement%d.txt", M_DIR, nm); 
    read_measurement(dim, M_FILE, mu, S, &T); 

    // time parameters
    clock_t start = clock();
    clock_t end;
    double rt;
    double dt = 0.001;

    // initialize cut
    double** Z = NULL; 
    double* W = NULL; 
    int N = 0; 
    if(CUT_M == 4){
        N = 2 * dim + (1 << dim)+ 1; 
        Z = init_m(N, dim); 
        W = init_v(N); 
        cut4(Z, W, mu, S, dim);
    }else if(CUT_M == 6){
        N = 2 * dim * dim + (1 << dim)+ 1; 
        Z = init_m(N, dim); 
        W = init_v(N); 
        cut6(Z, W, mu, S, dim);
    }else{
        fprintf(stderr, "Incorrect CUT order.");
        return -1; 
    }
    double a0 = 0.001; 
    double nu = dim * (dim + 1) * (dim + 2) / 6; 
    double z1_s = gsl_cdf_chisq_Qinv(a0, nu);
    double z2_s = gsl_cdf_ugaussian_Pinv(1.0 - (a0 / 2.0)); 
    double b1; 
    double b2; 
    calc_b1_b2(Z, W, N, CUT_M, dim, &b1, &b2);

    // file parameters
    char FILE_NAME[256];  
    char DIR_NAME[256];  
    FILE* file = NULL; 

    mkdir(R_DIR, 0777);
    snprintf(DIR_NAME, sizeof(DIR_NAME), "%s/sigmas/", R_DIR); 
    mkdir(DIR_NAME, 0777);
    snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/sigmas/%s_rv_cukf_%d.txt", R_DIR, M_MODEL, nm);
    file = fopen(FILE_NAME, "w");
    fprintf(file, "t, W, x [LU], y [LU], vx [LU/TU], vy [LU/TU], ..., b1, b2\n"); 
    fprintf(file, "%.3f ", 0.0); 
    for(int i = 0; i < N; i++){
        fprintf(file, "%.15e ", W[i]);
        for(int j = 0; j < dim; j++){
            fprintf(file, "%.15e ", Z[i][j]); 
        }
    }
    fprintf(file, "%.15f, %.15f\n", b1, b2);  
    printf("t: %.3f\n", 0.0);

    double t = 0.0;
    while((fabs(t - T) > 1e-4)&&(check_b1_b2(b1, b2, N, dim, z1_s, z2_s))){

        dt = fmin(dt, T - t); 

        for(int j = 0; j < N; j++){
            RK4(f, dim, Z[j], t, dt); 

            if(j == 0){
                fprintf(file, "%.3f ", t + dt); 
            }
            fprintf(file, "%.15e ", W[j]);
            for(int k = 0; k < dim; k++){
                fprintf(file, "%.15e ", Z[j][k]); 
            }
        }
        calc_b1_b2(Z, W, N, CUT_M, dim, &b1, &b2);
        fprintf(file, "%.15f, %.15f\n", b1, b2);  
        t += dt; 

        printf("t: %.3f, b1: %.15f, b2: %.15f\n", t, b1, b2);
    }
    end = clock();
    rt = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Total runtime: %.6f seconds\n", rt);

    // printing final 
    double* mu_cut = reconstruct_cut(Z, W, dim, 1, CUT_M);
    double* S_cut = reconstruct_cut(Z, W, dim, 2, CUT_M);
    double** S_mat = init_m(dim, dim);
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            S_mat[i][j] = S_cut[i * dim + j];
        }
    }
    print_vector(dim, mu_cut); 
    print_matrix(dim, dim, S_mat);
    free(mu_cut); free(S_cut);
    for(int i = 0; i < dim; i++){
        free(S_mat[i]); 
    }
    free(S_mat);
    printf("T: %.15f\n", T - t);

    // free memory
    fclose(file); 
    free(mu); 
    for(int i = 0; i < dim; i++){
        free(S[i]); 
    }
    free(S); 
    for(int i = 0; i < N; i++){
        free(Z[i]); 
    }
    free(Z); free(W); 
 
    return 0; 
}
