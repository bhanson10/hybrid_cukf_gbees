#include "../../src/utils.h"
#include "../../src/models.h"
#include "../../src/cukf.h"
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
    const char* HZ_s = "1-047";  
    double a0 = 0.003;
    int Npf = 5000; 
    double HZ_t = 1.047252575058097; // n = 5000, d = 6, alpha_0 = 0.003
    // double HZ_t = 1.041646662063038; // n = 10000, d = 6, alpha_0 = 0.003
    // double HZ_t = 1.029455469017074; // n = 50000, d = 6, alpha_0 = 0.003
    int CUT_M = 6; 
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
        snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s/cukfpf", orbit);
        mkdir(R_DIR, 0777);
        snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s/cukfpf/%s_%s_HZ_%s", orbit, len, pf_size, HZ_s);
        mkdir(R_DIR, 0777);
    }
    
    // time parameters
    double start = omp_get_wtime();
    double rt;  
    double tt = 0.0;

    // measurement update parameters
    double* y   = init_v(y_dim); 
    double** Sy = init_m(y_dim, y_dim);

    // initialize CUT
    double* Z = NULL; 
    double* W = NULL; 
    int N = 0; 
    if(CUT_M == 4){
        N = 2 * dim + (1 << dim)+ 1; 
        Z = init_v(N * dim); 
        W = init_v(N); 
        cut4(dim, Z, W, mu, S);
    }else if(CUT_M == 6){
        N = 2 * dim * dim + (1 << dim)+ 1; 
        Z = init_v(N * dim); 
        W = init_v(N); 
        cut6(dim, Z, W, mu, S);
    }else{
        fprintf(stderr, "Incorrect CUT order.");
        return -1; 
    } 
    double nu = dim * (dim + 1) * (dim + 2) / 6; 
    double z1_s = gsl_cdf_chisq_Qinv(a0, nu);
    double z2_s = gsl_cdf_ugaussian_Pinv(1.0 - (a0 / 2.0)); 
    double b1; 
    double b2; 
    bool cukf_b = 1; 
    int cukf_dir_c = 0;
    
    // PF parameters
    double HZ; 
    bool HZ_test = 0; 
    double* X; 
    double* Wpf = init_v(Npf); 
    bool switch_filter = 0; 
    int pf_dir_c = -1;
    double ess_t = 0.5 * Npf; 

    // Quadratic sampling parameters
    double eta = 1.0; 
    int HZ_step_min = 1; 
    int HZ_step_max = 50;
    int HZ_step = 1; 

    // file parameters
    char FILE_NAME[256];  
    char DIR_NAME[256];  
    FILE* file_cukf = NULL; 
    FILE* file_pf = NULL; 

    if(record){
        snprintf(DIR_NAME, sizeof(DIR_NAME), "%s/sigmas/", R_DIR); 
        mkdir(DIR_NAME, 0777);
        snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/sigmas/%s_rv_cukf_0.txt", R_DIR, orbit);
        file_cukf = fopen(FILE_NAME, "w");
        fprintf(file_cukf, "t, W, x [LU], y [LU], z [LU], vx [LU/TU], vy [LU/TU], vz [LU/TU], ...\n"); 
        fprintf(file_cukf, "%.3f ", 0.0); 
        for(int i = 0; i < N; i++){
            fprintf(file_cukf, "%.15e ", W[i]);
            for(int j = 0; j < dim; j++){
                fprintf(file_cukf, "%.15e ", Z[i * dim + j]); 
            }
        } 
        fprintf(file_cukf, "\n"); 
    }
    calc_bt1_bt2(N, dim, CUT_M, Z, W, mu, S, &b1, &b2);
    printf("t: %.3f, b1: %.12f, b2: %.12f, PF mode: %d\n", 0.0, b1, b2, (int) !cukf_b);

    int nm = count_txt_files(M_DIR);
    for(int i = 1; i <= nm; i++){

        // Integrate until we reach tend or an error occurs
        double t = 0.0; 
        double dt = 0.001;
        int step_count = 1; 
        while(fabs(t - T) > 1e-4){

            dt = fmin(dt, T - t); 
            
            if(cukf_b){
                for(int j = 0; j < N; j++){
                    RK4(f, dim, &Z[j * dim], coef_f, t, dt); 

                    if(record){
                        if(j == 0){
                            fprintf(file_cukf, "%.3f ", tt + t + dt); 
                        }
                        fprintf(file_cukf, "%.15e ", W[j]);
                        for(int k = 0; k < dim; k++){
                            fprintf(file_cukf, "%.15e ", Z[j * dim + k]); 
                        }
                    }
                }
                if(record) fprintf(file_cukf, "\n"); 
            }else{

                if(!HZ_test){
                    for(int j = 0; j < N; j++){
                        RK4(f, dim, &Z[j * dim], coef_f, t, dt); 
                        
                        if(record){
                            if(j == 0){
                                fprintf(file_cukf, "%.3f ", tt + t + dt); 
                            }
                            fprintf(file_cukf, "%.15e ", W[j]);
                            for(int k = 0; k < dim; k++){
                                fprintf(file_cukf, "%.15e ", Z[j * dim + k]); 
                            }
                        }
                    }
                    if(record) fprintf(file_cukf, "\n");  
                }

                for(int j = 0; j < Npf; j++){
                    RK4(f, dim, &X[j * dim], coef_f, t, dt); 

                    if(record){
                        if(j == 0){
                            fprintf(file_pf, "%.3f ", tt + t + dt); 
                        }
                        fprintf(file_pf, "%.15e ", Wpf[j]);
                        for(int k = 0; k < dim; k++){
                            fprintf(file_pf, "%.15e ", X[j * dim + k]); 
                        }
                    }
                }
                if(record) fprintf(file_pf, "\n");  
            }
            t += dt; 
            
            if(!HZ_test){
                calc_b1_b2(N, dim, CUT_M, Z, W, mu, S, &b1, &b2);
                if(cukf_b){
                    if(!check_b1_b2(b1, b2, N, dim, z1_s, z2_s)) switch_filter = 1; 
                }else{
                    if(check_b1_b2(b1, b2, N, dim, z1_s, z2_s)) switch_filter = 1; 
                }
                printf("t: %.3f, b1: %.12f, b2: %.12f, PF mode: %d\n", tt + t, b1, b2, (int) !cukf_b);
            }else{
                if((step_count % HZ_step == 0)||(fabs(t - T) < 1e-4)){
                    HZ = calc_HZ_fast(Npf, dim, X, NULL, mu, S);
                    if(HZ < HZ_t) switch_filter = 1; 
                    HZ_step = fmax(HZ_step_min, (int)ceil(HZ_step_max / (1 + eta / fabs(HZ - HZ_t))));
                    step_count = 0;
                    printf("t: %.3f, HZ: %.12f, PF mode: %d\n", tt + t, HZ, (int) !cukf_b);
                }else{
                    printf("t: %.3f, HZ: NaN, PF mode: %d\n", tt + t, (int) !cukf_b);
                }
                step_count++;
            }
        
            // switching
            if(switch_filter){
                if(cukf_b){
                    cukf_b = 0; switch_filter = 0; 

                    // sampling from last CUKF
                    X = mvnrnd(Npf, dim, mu, S);
                    for(int j = 0; j < Npf; j++){
                        Wpf[j] = 1.0 / Npf; 
                    }
                    
                    if(record){
                        snprintf(DIR_NAME, sizeof(DIR_NAME), "%s/particles/", R_DIR); 
                        if(pf_dir_c == -1) mkdir(DIR_NAME, 0777);
                        snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/particles/%s_rv_pf_%d.txt", R_DIR, orbit, ++pf_dir_c);
                        file_pf = fopen(FILE_NAME, "w");
                        fprintf(file_pf, "t, W, x [LU], y [LU], z [LU], vx [LU/TU], vy [LU/TU], vz [LU/TU], ...\n"); 
                        fprintf(file_pf, "%.3f ", tt + t); 
                        for(int j = 0; j < Npf; j++){
                            fprintf(file_pf, "%.15e ", Wpf[j]);
                            for(int k = 0; k < dim; k++){
                                fprintf(file_pf, "%.15e ", X[j * dim + k]); 
                            }
                        }
                        fprintf(file_pf, "\n");  
                    }

                }else{
                    cukf_b = 1; switch_filter = 0; 

                    if(HZ_test){
                        HZ_test = 0; 

                        fclose(file_pf); file_pf = NULL; 
                        fclose(file_cukf); file_cukf = NULL; 

                        // initializing CUKF from last pf 
                        if(CUT_M == 4){
                            cut4(dim, Z, W, mu, S);
                        }else if(CUT_M == 6){
                            cut6(dim, Z, W, mu, S);
                        }
                        
                        if(record){
                            snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/sigmas/%s_rv_cukf_%d.txt", R_DIR, orbit, ++cukf_dir_c);
                            file_cukf = fopen(FILE_NAME, "w");
                            fprintf(file_cukf, "t, W, x [LU], y [LU], z [LU], vx [LU/TU], vy [LU/TU], vz [LU/TU], ...\n"); 
                            fprintf(file_cukf, "%.3f ", tt + t); 
                            for(int j = 0; j < N; j++){
                                fprintf(file_cukf, "%.15e ", W[j]);
                                for(int k = 0; k < dim; k++){
                                    fprintf(file_cukf, "%.15e ", Z[j * dim + k]); 
                                }
                            }
                            fprintf(file_cukf, "\n");  
                        }
                        calc_b1_b2(N, dim, CUT_M, Z, W, mu, S, &b1, &b2);
                        printf("t: %.3f, b1: %.12f, b2: %.12f, PF mode: %d\n", tt + t + dt, b1, b2, (int) !cukf_b);
                    }
                }
            }
        }
        tt += T;

        // measurement update
        if(i != nm){

            snprintf(M_FILE, sizeof(M_FILE), "%s/measurement%d.txt", M_DIR, i);
            read_measurement(y_dim, M_FILE, y, Sy, &T); 

            if(cukf_b){
                measurement_update_cukf(h, &Z, &W, y, coef_h, Sy, tt, N, dim, y_dim, CUT_M);

                fclose(file_cukf); file_cukf = NULL; 
                
                if(record){
                    snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/sigmas/%s_rv_cukf_%d.txt", R_DIR, orbit, ++cukf_dir_c);
                    file_cukf = fopen(FILE_NAME, "w");
                    fprintf(file_cukf, "t, W, x [LU], y [LU], z [LU], vx [LU/TU], vy [LU/TU], vz [LU/TU], ...\n"); 
                    fprintf(file_cukf, "%.3f ", tt); 
                    for(int j = 0; j < N; j++){
                        fprintf(file_cukf, "%.15e ", W[j]);
                        for(int k = 0; k < dim; k++){
                            fprintf(file_cukf, "%.15e ", Z[j * dim + k]); 
                        }
                    }
                    fprintf(file_cukf, "\n");  
                }
                calc_b1_b2(N, dim, CUT_M, Z, W, mu, S, &b1, &b2);
                printf("t: %.3f, b1: %.12f, b2: %.12f, PF mode: %d\n", tt + t + dt, b1, b2, (int) !cukf_b);

            }else{

                measurement_update_pf(h, X, Wpf, y, coef_h, Sy, tt, ess_t, Npf, dim, y_dim);
                HZ = calc_HZ_fast(Npf, dim, X, NULL, mu, S);

                if(HZ < HZ_t){
                    cukf_b = 1; switch_filter = 0; 

                    fclose(file_pf); file_pf = NULL; 
                    fclose(file_cukf); file_cukf = NULL; 

                    // generate CUKF from new mu and S
                    if(CUT_M == 4){
                        cut4(dim, Z, W, mu, S);
                    }else if(CUT_M == 6){
                        cut6(dim, Z, W, mu, S);
                    }

                    if(record){
                        snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/sigmas/%s_rv_cukf_%d.txt", R_DIR, orbit, ++cukf_dir_c);
                        file_cukf = fopen(FILE_NAME, "w");
                        fprintf(file_cukf, "t, W, x [LU], y [LU], z [LU], vx [LU/TU], vy [LU/TU], vz [LU/TU], ...\n"); 
                        fprintf(file_cukf, "%.3f ", tt + t); 
                        for(int j = 0; j < N; j++){
                            fprintf(file_cukf, "%.15e ", W[j]);
                            for(int k = 0; k < dim; k++){
                                fprintf(file_cukf, "%.15e ", Z[j * dim + k]); 
                            }
                        }
                        fprintf(file_cukf, "\n"); 
                    }
                    calc_b1_b2(N, dim, CUT_M, Z, W, mu, S, &b1, &b2);
                    printf("t: %.3f, b1: %.12f, b2: %.12f, PF mode: %d\n", tt + t + dt, b1, b2, (int) !cukf_b);

                }else{
                    HZ_test = 1; switch_filter = 0; 
                    for(int j = 0; j < Npf; j++){
                        Wpf[j] = 1.0 / Npf; 
                    }

                    fclose(file_pf); file_pf = NULL; 
                    fclose(file_cukf); file_cukf = NULL; 
                    
                    if(record){
                        snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/particles/%s_rv_pf_%d.txt", R_DIR, orbit, ++pf_dir_c);
                        file_pf = fopen(FILE_NAME, "w");
                        fprintf(file_pf, "t, W, x [LU], y [LU], z [LU], vx [LU/TU], vy [LU/TU], vz [LU/TU], ...\n"); 
                        fprintf(file_pf, "%.3f ", tt); 
                        for(int j = 0; j < Npf; j++){
                            fprintf(file_pf, "%.15e ", Wpf[j]);
                            for(int k = 0; k < dim; k++){
                                fprintf(file_pf, "%.15e ", X[j * dim + k]); 
                            }
                        }
                        fprintf(file_pf, "\n");  
                    }
                    printf("t: %.3f, HZ: %.12f, PF mode: %d\n", tt, HZ, (int) !cukf_b);
                }
                
            }
        }
    }
    double end = omp_get_wtime();
    rt = end - start;
    printf("Total runtime: %.6f seconds\n", rt);

    // Clean up memory
    if (file_pf != NULL){
        fclose(file_pf);
    }
    if (file_cukf != NULL){
        fclose(file_pf);
    }
    free(mu); 
    for(int i = 0; i < dim; i++){
        free(S[i]); 
    }
    free(S); 
    free(Z); free(W); 
    free(X); free(Wpf); 

    printf("Complete.\n");

    return 0; 
}
