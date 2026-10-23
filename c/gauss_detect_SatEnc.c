#include "./src/utils.h"
#include "./src/models.h"
#include "./src/ukf.h"
#include "./src/cukf.h"

int main(int argc, char* argv[]){

    // inputs
    fn f = CR3BP;
    int dim = 6;    
    int N = 50; 
    int Nfam = 10; 
    char fams[10][100] = {"DPO", "DRAGONFLY_NORTHERN", "DRAGONFLY_SOUTHERN", "DRO", "L1_HALO_NORTHERN", "L1_HALO_SOUTHERN", "L2_BUTTERFLY_NORTHERN", "L2_BUTTERFLY_SOUTHERN", "L2_HALO_NORTHERN", "L2_HALO_SOUTHERN"};
    double sig_r = 1;       // km
    double sig_v = 0.01e-3; // km/s
    double alpha_0 = 0.003;
    double NED_t = 0.1296; 
    double coef[1] = {1.901109735892602E-7}; // Saturn-Enceladus system

    // parameters
    double LU = 238529; 
    double TU = 18913;

    // MC variables
    int Nmc = 5000;  // number of samples
    double b_HZ = (1.0 / sqrt(2.0)) * pow((2.0 * dim + 1.0) / 4.0, 1.0 / (dim + 4.0)) * pow(Nmc, 1.0 / (dim + 4.0));
    double a_HZ = 1.0 + 2.0 * pow(b_HZ, 2.0);
    double wb = (1.0 + pow(b_HZ, 2.0)) * (1.0 + 3.0 * pow(b_HZ, 2.0));
    double mu_h = 1.0 - pow(a_HZ, -dim / 2.0) * (1.0 + dim * pow(b_HZ, 2.0) / a_HZ + (dim * (dim + 2.0) * pow(b_HZ, 4.0)) / (2.0 * pow(a_HZ, 2.0)));
    double s2_h = 2.0 * pow(1.0 + 4.0 * pow(b_HZ, 2.0), -dim / 2.0)
                + 2.0 * pow(a_HZ, -dim) * (1.0 + (2.0 * dim * pow(b_HZ, 4.0)) / pow(a_HZ, 2.0) + (3.0 * dim * (dim + 2.0) * pow(b_HZ, 8.0)) / (4.0 * pow(a_HZ, 4.0)))
                - 4.0 * pow(wb, -dim / 2.0) * (1.0 + (3.0 * dim * pow(b_HZ, 4.0)) / (2.0 * wb) + (dim * (dim + 2.0) * pow(b_HZ, 8.0)) / (2.0 * pow(wb, 2.0)));
    double ratio = s2_h / pow(mu_h, 2.0);
    double norminv_val = gsl_cdf_ugaussian_Pinv(1.0 - alpha_0);
    double HZ_t = mu_h * pow(1.0 + ratio, -0.5) * exp(norminv_val * sqrt(log(1.0 + ratio)));

    printf("HZ_t: %.4f\n", HZ_t); 

    // UKF variables
    int Nukf = 2 * dim + 1; // number of weighted sigma points
    double a = 0.8;         // alpha
    double b = 2.0;         // beta
    double k = 0.0;         // kappa
    
    // CUKF6 variables
    int CUT_M = 6; 
    int Ncukf = 2 * dim * dim + (1 << dim) + 1; 
    double nu = dim * (dim + 1) * (dim + 2) / 6; 
    double z1_t = gsl_cdf_chisq_Qinv(alpha_0, nu);
    double z2_t = gsl_cdf_ugaussian_Pinv(1.0 - (alpha_0 / 2.0)); 
    double b1_t = z1_t / ((dim + 1.0) * (Ncukf + 1.0) * (Ncukf + 3.0) / (6.0 * ((Ncukf + 1.0) * (dim + 1.0) - 6.0)));
    double b2_t = z2_t * sqrt(8.0 * dim * (dim + 2.0) / Ncukf) + dim * (dim + 2.0);

    printf("z1_t: %.4f\n", z1_t); 
    printf("z2_t: %.4f\n", z2_t); 

    return 0; 

    printf("GVT (TU) order: HZ, b, NED, bt\n");

    #pragma omp parallel for collapse(2) schedule(dynamic, 1)
    for(int i = 0; i < Nfam; i++){
        for(int j = 0; j < N; j++){

            // Thread-private loop variables
            double HZ = 0.0, b1 = 0.0, b2 = 0.0, NED = 0.0, bt1 = 0.0, bt2 = 0.0; 
            double t_HZ = 0.0, t_b = 0.0, t_NED = 0.0, t_bt = 0.0; 

            // Thread-private state pointers
            double* mu    = init_v(dim);
            double* X     = NULL;
            double* Zukf  = init_v(Nukf * dim);
            double* Wukf  = init_v(Nukf * 2);
            double* Zcukf = init_v(Ncukf * dim); 
            double* Wcukf = init_v(Ncukf); 
            
            bool HZ_b = true;
            bool NED_b = true;
            bool b_b = true;
            bool bt_b = true;   
            double T = 0; 

            char FILE_NAME[256];
            snprintf(FILE_NAME, sizeof(FILE_NAME), "/Users/behanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/hybrid_cr3bp/results/python/families/trajs/%s/SatEnc_CR3BP_%d_rv.txt", fams[i], j); 
            read_measurement_2_p(dim, FILE_NAME, mu, &T); 

            double** S = init_m(dim, dim);
            for(int k = 0; k < dim; k++){
                for(int l = 0; l < dim; l++){
                    if(k == l && k < 3) S[k][l] = pow(sig_r / LU, 2.0); 
                    else if(k == l && k >= 3) S[k][l] = pow(sig_v / (LU / TU), 2.0); 
                    else S[k][l] = 0.0; 
                }
            }

            unsigned int seed = (unsigned int)(1469598103U ^ ((unsigned int)i * 73856093U) ^ ((unsigned int)j * 19349663U));
            X = mvnrnd_p(Nmc, dim, mu, S, &seed);
            ut(dim, Zukf, Wukf, mu, S, a, b, k);
            cut6(dim, Zcukf, Wcukf, mu, S);

            // snprintf(FILE_NAME, sizeof(FILE_NAME), "../results/gauss_detect/SatEnc/%s/SatEnc_CR3BP_gd_%d.txt", fams[i], j); 
            // FILE* file = fopen(FILE_NAME, "w");
            // fprintf(file, "t (TU), t (hr), t (T), HZ, b1, b2, NED, bt1, bt2\n"); 

            double t = 0.0;
            double dt = 0.001;
            while(HZ_b || b_b || NED_b || bt_b){
            // while (fabs(t - T) > 1e-4){

                if(HZ_b){
                    HZ = calc_HZ_fast(Nmc, dim, X, NULL, mu, S);
                    if(HZ > HZ_t){
                        HZ_b = false;
                        t_HZ = t; 
                    }
                }

                if(b_b){
                    b1 = calc_b1(Nmc, dim, X, mu, S);
                    if(b1 > b1_t){
                        b_b = false;
                        t_b = t; 
                    }else{
                        b2 = calc_b2(Nmc, dim, X, mu, S);
                        if((b2 > b2_t) || (b2 < -b2_t)){
                            b_b = false;
                            t_b = t; 
                        }
                    }
                }

                if(NED_b){
                    NED = calc_NED(dim, Zukf, Wukf, mu, S);
                    if(NED > NED_t){
                        NED_b = false;
                        t_NED = t;
                    }
                }

                if(bt_b){
                    bt1 = calc_bt1(Ncukf, dim, CUT_M, Zcukf, Wcukf, mu, S);
                    if(bt1 > b1_t){
                        bt_b = false;
                        t_bt = t; 
                    }else{
                        bt2 = calc_bt2(Ncukf, dim, CUT_M, Zcukf, Wcukf, mu, S);
                        if((bt2 > b2_t) || (bt2 < -b2_t)){
                            bt_b = false;
                            t_bt = t; 
                        }
                    }
                }


                if(!HZ_b && !b_b && !NED_b && !bt_b) break;

                // propagate particles
                // dt = fmin(dt, T - t); 
                for(int k = 0; k < Nmc; k++){
                    RK4(f, dim, &X[k * dim], coef, t, dt); 
                }
                for(int k = 0; k < Nukf; k++){
                    RK4(f, dim, &Zukf[k * dim], coef, t, dt); 
                }
                for(int k = 0; k < Ncukf; k++){
                    RK4(f, dim, &Zcukf[k * dim], coef, t, dt); 
                }
                t += dt; 
            }

            #pragma omp critical
            {
                printf("Family: %s, IC: %d, Period (TU): %.10f - %.3f %.3f %.3f %.3f\n", fams[i], j, T, t_HZ, t_b, t_NED, t_bt);
            }

            // Clean up memory
            free(mu); mu = NULL;
            for(int k = 0; k < dim; k++){
                free(S[k]); 
            }
            free(S); 
            free(X); X = NULL; 
            free(Zukf); free(Wukf); Zukf = NULL; Wukf = NULL;
            free(Zcukf); free(Wcukf); Zcukf = NULL; Wcukf = NULL;
            // fclose(file); file = NULL;
        }
    }

    return 0; 
}
