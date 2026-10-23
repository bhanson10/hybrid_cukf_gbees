#include "../../src/utils.h"
#include "../../src/models.h"
#include "../../src/cukf.h"
#include "../../src/pme.h"

int main(int argc, char* argv[]){

    // inputs
    fn f = PCR3BP;
    int dim = 4;              
    fn h = LOS_2;
    int y_dim = 4;
    const char* M_MODEL = "PCR3BP";
    int CUT_M = 6;
    bool pme = false; 
    int PME_M = 3; 
    int n_quad = 50; 
    int n_eval = 75;
    
    // parameters
    double* mu  = init_v(dim); 
    double** S  = init_m(dim, dim);
    double T    = 0; 
    
    // directories
    char M_DIR[256];  
    snprintf(M_DIR, sizeof(M_DIR), "../../measurements/%s", M_MODEL);
    char R_DIR[256];  
    snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s", M_MODEL);
    mkdir(R_DIR, 0777);
    if(pme){
        snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s/cukf_pme", M_MODEL);
    }else{
        snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s/cukf", M_MODEL);
    }
    mkdir(R_DIR, 0777);
    char M_FILE[256];  
    snprintf(M_FILE, sizeof(M_FILE), "%s/measurement0.txt", M_DIR); 
    read_measurement(dim, M_FILE, mu, S, &T); 

    // time parameters
    clock_t start = clock();
    clock_t end;
    double rt;
    double tt = 0.0;

    // measurement update parameters
    double* y   = init_v(y_dim); 
    double** Sy = init_m(y_dim, y_dim);

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

    // PME reconstruction
    double* X = NULL; 
    double* P = NULL; 
    int num_pts = 1;
    for(int i = 0; i < dim; i++) num_pts *= n_eval;

    // file parameters
    char FILE_NAME[256];  
    char DIR_NAME[256];  
    FILE* file = NULL; 

    snprintf(DIR_NAME, sizeof(DIR_NAME), "%s/sigmas/", R_DIR); 
    mkdir(DIR_NAME, 0777);
    if(pme){
        snprintf(DIR_NAME, sizeof(DIR_NAME), "%s/pdfs/", R_DIR); 
        mkdir(DIR_NAME, 0777);
        snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/sigmas/%s_rv_cukf_pme_0.txt", R_DIR, M_MODEL);
    }else{
        snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/sigmas/%s_rv_cukf_0.txt", R_DIR, M_MODEL);
    }
    file = fopen(FILE_NAME, "w");
    fprintf(file, "t, W, x [LU], y [LU], vx [LU/TU], vy [LU/TU], ...\n"); 
    fprintf(file, "%.6f ", 0.0); 
    for(int i = 0; i < N; i++){
        fprintf(file, "%.15e ", W[i]);
        for(int j = 0; j < dim; j++){
            fprintf(file, "%.15e ", Z[i][j]); 
        }
    }
    fprintf(file, "\n");  
    printf("t: %.6f\n", 0.0);
    
    int nm = count_txt_files(M_DIR);
    for(int i = 1; i <= nm; i++){

        double t = 0.0;
        double dt = 0.001;
        while (fabs(t - T) > 1e-4){

            dt = fmin(dt, T - t); 

            for(int j = 0; j < N; j++){
                RK4(f, dim, Z[j], tt + t, dt); 

                if(j == 0){
                    fprintf(file, "%.6f ", tt + t + dt); 
                }
                fprintf(file, "%.15e ", W[j]);
                for(int k = 0; k < dim; k++){
                    fprintf(file, "%.15e ", Z[j][k]); 
                }
            }
            fprintf(file, "\n");  
            t += dt; 

            printf("t: %.6f\n", tt + t);
        }
        tt += T;

        if(pme){
            pme_cut2pdf(Z, W, N, CUT_M, dim, PME_M, n_quad, n_eval, &X, &P);
            snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/pdfs/%s_rv_pdf_%d.bin", R_DIR, M_MODEL, i - 1);
            write_pdf(FILE_NAME, t, num_pts, dim, X, P);
        }

        // measurement update
        if(i != nm){

            snprintf(M_FILE, sizeof(M_FILE), "%s/measurement%d.txt", M_DIR, i);
            read_measurement(y_dim, M_FILE, y, Sy, &T);  

            if(pme){
                measurement_update_cukf_pme(h, &Z, &W, X, P, num_pts, y, Sy, tt, N, dim, y_dim, CUT_M);
            }else{
                measurement_update_cukf(h, &Z, &W, y, Sy, tt, N, dim, y_dim, CUT_M);
            }
            
            fclose(file); file = NULL; 
            if(pme){
                snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/sigmas/%s_rv_cukf_pme_%d.txt", R_DIR, M_MODEL, i);
            }else{
                snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/sigmas/%s_rv_cukf_%d.txt", R_DIR, M_MODEL, i);
            }
            file = fopen(FILE_NAME, "w");
            fprintf(file, "t, W, x [LU], y [LU], vx [LU/TU], vy [LU/TU], ...\n"); 
            fprintf(file, "%.6f ", tt); 
            for(int j = 0; j < N; j++){
                fprintf(file, "%.15e ", W[j]);
                for(int k = 0; k < dim; k++){
                    fprintf(file, "%.15e ", Z[j][k]); 
                }
            }   
            fprintf(file, "\n");  
            printf("t: %.6f\n", tt);
        }
    }
    end = clock();
    rt = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Total runtime: %.6f seconds\n", rt);

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
