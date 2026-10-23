#include "../../src/utils.h"
#include "../../src/models.h"
#include "../../src/pf.h"

int main(int argc, char* argv[]){

    // inputs
    fn f = CR3BP;
    int dim = 6;    
    fn h = R_RR; 
    int y_dim = 2;     
    const char* orbit = "NRHO"; 
    const char* len = "30d"; 
    const char* pf_size = "5000";
    int N = 5000; 
    bool record = false; 
    double coef_f[1] = {0.012150585609624}; // Earth-Moon system
    double coef_h[1] = {0.83691513};        // Earth-Moon L1 x-value

    // parameters
    srand(time(NULL)); // Seed the RNG with current time            
    double* mu  = init_v(dim); 
    double** S  = init_m(dim, dim);
    double T    = 0; 
    char M_DIR[256];
    snprintf(M_DIR, sizeof(M_DIR), "../../measurements/%s/cukfpf/%s", orbit, len);
    char M_FILE[256];  
    snprintf(M_FILE, sizeof(M_FILE), "%s/measurement0.txt", M_DIR); 
    read_measurement(dim, M_FILE, mu, S, &T);  

    // directories
    char R_DIR[256];  
    if(record){
        snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s", orbit);
        mkdir(R_DIR, 0777);
        snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s/pf", orbit);
        mkdir(R_DIR, 0777);
        snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s/pf/%s_%s", orbit, len, pf_size);
        mkdir(R_DIR, 0777);
    }
    
    // time parameters
    double start = omp_get_wtime();
    double rt;
    double tt = 0.0;

    // measurement update parameters
    double* y   = init_v(y_dim); 
    double** Sy = init_m(y_dim, y_dim);

    // PF parameters
    double* X = mvnrnd(N, dim, mu, S);
    double* W = init_v(N); 
    double ess_t = 0.5 * N; 

    // file parameters
    char FILE_NAME[256];  
    char DIR_NAME[256];  
    FILE* file = NULL; 

    if(record){
        snprintf(DIR_NAME, sizeof(DIR_NAME), "%s/particles/", R_DIR); 
        mkdir(DIR_NAME, 0777);
        snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/particles/%s_rv_pf_0.txt", R_DIR, orbit);
        file = fopen(FILE_NAME, "w");
        fprintf(file, "t, W, x [LU], y [LU], z [LU], vx [LU/TU], vy [LU/TU], vz [LU/TU], ...\n"); 
        fprintf(file, "%.3f ", 0.0); 
        for(int i = 0; i < N; i++){
            W[i] = 1.0 / N; 
            fprintf(file, "%.15e ", W[i]);
            for(int j = 0; j < dim; j++){
                fprintf(file, "%.15e ", X[i * dim + j]); 
            }
        }
        fprintf(file, "\n");  
    }
    printf("t: %.3f\n", 0.0);

    int nm = count_txt_files(M_DIR);
    for(int i = 1; i <= nm; i++){

        double t = 0.0; 
        double dt = 0.001;
        while (fabs(t - T) > 1e-4){

            dt = fmin(dt, T - t); 

            for(int j = 0; j < N; j++){
                RK4(f, dim, &X[j * dim], coef_f, t, dt); 

                if(record){
                    if(j == 0){
                        fprintf(file, "%.3f ", tt + t + dt); 
                    }
                    fprintf(file, "%.15e ", W[j]);
                    for(int k = 0; k < dim; k++){
                        fprintf(file, "%.15e ", X[j * dim + k]); 
                    }
                }
            }
            if(record) fprintf(file, "\n");  
            t += dt; 

            printf("t: %.3f\n", tt + t);
        }
        tt += T;

        // measurement update
        if(i != nm){

            snprintf(M_FILE, sizeof(M_FILE), "%s/measurement%d.txt", M_DIR, i);
            read_measurement(y_dim, M_FILE, y, Sy, &T); 

            measurement_update_pf(h, X, W, y, coef_h, Sy, tt, ess_t, N, dim, y_dim);
            mu = mean(N, dim, X, NULL); 
            S = cov(N, dim, X, NULL);

            fclose(file); file = NULL; 
            if(record){
                snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/particles/%s_rv_pf_%d.txt", R_DIR, orbit, i);
                file = fopen(FILE_NAME, "w");
                fprintf(file, "t, W, x [LU], y [LU], z [LU], vx [LU/TU], vy [LU/TU], vz [LU/TU], ...\n"); 
                fprintf(file, "%.3f ", tt); 
                for(int j = 0; j < N; j++){
                    fprintf(file, "%.15e ", W[j]);
                    for(int k = 0; k < dim; k++){
                        fprintf(file, "%.15e ", X[j * dim + k]); 
                    }
                }   
                fprintf(file, "\n");  
            }
            printf("t: %.3f\n", tt);
        }
    }
    double end = omp_get_wtime();
    rt = end - start;
    printf("Total runtime: %.6f seconds\n", rt);

    // Clean up memory
    fclose(file); 
    free(mu); 
    for(int i = 0; i < dim; i++){
        free(S[i]); 
    }
    free(S); 
    free(X); free(W); 

    printf("Complete.\n");

    return 0; 
}
