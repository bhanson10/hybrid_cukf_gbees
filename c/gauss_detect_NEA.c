#include "rebound.h"
#include "assist.h"
#include "./src/utils.h"
#include "./src/utils_assist.h"
#include "./src/ukf.h"
#include "./src/cukf.h"

int main(int argc, char* argv[]){

    // inputs
    int dim = 6;    
    int N   = 1000; 
    double alpha_0 = 0.003;
    double NED_t = 0.14668; 
    
    // assist variables
    struct assist_ephem* ephem = assist_ephem_create(
            "/Users/behanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/rebound-assist/data/linux_p1550p2650.440",
            "/Users/behanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/rebound-assist/data/sb441-n16.bsp");
    if (!ephem){
        printf("Cannot create ephemeris structure.\n");
        exit(-1);
    }

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

    // UKF variables
    int Nukf = 2 * dim + 1;  // number of weighted sigma points
    double a = 0.8;          // alpha
    double b = 2.0;          // beta
    double k = 0.0;          // kappa
    
    // CUKF6 variables
    int CUT_M = 6; 
    int Ncukf = 2 * dim * dim + (1 << dim)+ 1; 
    double nu = dim * (dim + 1) * (dim + 2) / 6; 
    double z1_t = gsl_cdf_chisq_Qinv(alpha_0, nu);
    double z2_t = gsl_cdf_ugaussian_Pinv(1.0 - (alpha_0 / 2.0)); 
    double b1_t = z1_t / ((dim + 1.0) * (Ncukf + 1.0) * (Ncukf + 3.0) / (6.0 * ((Ncukf + 1.0) * (dim + 1.0) - 6.0)));
    double b2_t = z2_t * sqrt(8.0 * dim * (dim + 2.0) / Ncukf) + dim * (dim + 2.0);

    // misc. variables
    double t_thresh = 1; 
    double max_dt = 50; // Maximum step size of 50 days  

    double start = omp_get_wtime();

    printf("GVT (d) order: HZ, b, NED, bt\n");

    #pragma omp parallel for schedule(dynamic)
    for(int i = 201; i <= N; i++){

        // Thread-private loop variables
        double HZ = 0.0, b1 = 0.0, b2 = 0.0, NED = 0.0, bt1 = 0.0, bt2 = 0.0; 
        double t_HZ = 0.0, t_b = 0.0, t_NED = 0.0, t_bt = 0.0; 
        
        // Thread-private state pointers
        double* mu    = init_v(dim);
        double** S    = init_m(dim, dim);
        double* X     = NULL;
        double* Zukf  = init_v(Nukf * dim);
        double* Wukf  = init_v(Nukf * 2);
        double* Zcukf = init_v(Ncukf * dim); 
        double* Wcukf = init_v(Ncukf); 

        bool HZ_b = true;
        bool NED_b = true;
        bool b_b = true;
        bool bt_b = true;   
        double JD_0 = 0;     
        double T    = 0; 

        struct reb_simulation* r = NULL;
        struct assist_extras* ax = NULL;

        char FILE_NAME[256];
        snprintf(FILE_NAME, sizeof(FILE_NAME), "/Users/behanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Conferences/ISSFD_2026/ISSFD2026/c/NEA/meas/meas/NEA_rv_IC_%d.txt", i); 

        FILE* file_exists = fopen(FILE_NAME, "r");
        if(file_exists){
            read_measurement_p(dim, FILE_NAME, mu, S, &JD_0, &T); 
            T /= 5.0; // periods were pre-multiplied by 5 in the previous file
            double A = frob_norm(dim, dim, S); 

            unsigned int seed = (unsigned int)(1469598103U ^ ((unsigned int)i * 73856093U) ^ ((unsigned int)i * 19349663U));
            X = mvnrnd_p(Nmc, dim, mu, S, &seed);
            ut(dim, Zukf, Wukf, mu, S, a, b, k);
            cut6(dim, Zcukf, Wcukf, mu, S);

            // create simulation
            r = reb_simulation_create(); // Create a REBOUND simulation
            r->ri_ias15.adaptive_mode = 1;                 
            r->ri_ias15.min_dt = 300.0 / 86400; // Minimum step size of 5 minutes
            ax = assist_attach(r, ephem);       // Attach the ASSIST framework
            ax->gr_eih_sources = 11;            // Turn on general relativistic corrections for all bodies
            r->t = JD_0 - ephem->jd_ref;        // Initial time, relative to ephem->jd_ref
            
            // add NEA MC samples
            for(int j = 0; j < Nmc; j++){
                reb_simulation_add_fmt(r, "x y z vx vy vz",
                    X[j * dim + 0],  // x in AU
                    X[j * dim + 1],  // y
                    X[j * dim + 2],  // z
                    X[j * dim + 3],  // vx in AU/day
                    X[j * dim + 4],  // vy
                    X[j * dim + 5]); // vz
            }

            // add NEA UKF sigmas
            for(int j = 0; j < Nukf; j++){
                reb_simulation_add_fmt(r, "x y z vx vy vz",
                    Zukf[j * dim + 0],  // x in AU
                    Zukf[j * dim + 1],  // y
                    Zukf[j * dim + 2],  // z
                    Zukf[j * dim + 3],  // vx in AU/day
                    Zukf[j * dim + 4],  // vy
                    Zukf[j * dim + 5]); // vz
            }

            // add NEA CUKF sigmas
            for(int j = 0; j < Ncukf; j++){
                reb_simulation_add_fmt(r, "x y z vx vy vz",
                    Zcukf[j * dim + 0],  // x in AU
                    Zcukf[j * dim + 1],  // y
                    Zcukf[j * dim + 2],  // z
                    Zcukf[j * dim + 3],  // vx in AU/day
                    Zcukf[j * dim + 4],  // vy
                    Zcukf[j * dim + 5]); // vz
            }

            double t_last = r->t - t_thresh - 1.0;
            double t_end = r->t + 50 * T; // maximum propagation period of 50 periods
            while((HZ_b || b_b || NED_b || bt_b)&&(r->t <= t_end && r->status <= 0)){
            // while(r->t <= t_end && r->status <= 0){

                if((r->t - t_last) > t_thresh){
                    if(HZ_b){
                        HZ = calc_HZ_fast_particles(Nmc, dim, r->particles, NULL, mu, S);
                        if(HZ > HZ_t){
                            HZ_b = false;
                            t_HZ = r->t - (JD_0 - ephem->jd_ref); 
                        }
                    }

                    if(b_b){
                        b1 = calc_b1_particles(Nmc, dim, r->particles, mu, S);
                        if(b1 > b1_t){
                            b_b = false;
                            t_b = r->t - (JD_0 - ephem->jd_ref); 
                        }else{
                            b2 = calc_b2_particles(Nmc, dim, r->particles, mu, S);
                            if((b2 > b2_t) || (b2 < -b2_t)){
                                b_b = false;
                                t_b = r->t - (JD_0 - ephem->jd_ref); 
                            }
                        }
                    }

                    if(NED_b){
                        NED = calc_NED_particles(dim, &r->particles[Nmc], Wukf, mu, S);
                        if(NED > NED_t){
                            NED_b = false;
                            t_NED = r->t - (JD_0 - ephem->jd_ref);
                        }
                    }
                    
                    if(bt_b){
                        bt1 = calc_bt1_particles(Ncukf, dim, CUT_M, &r->particles[Nmc + Nukf], Wcukf, mu, S);
                        if(bt1 > b1_t){
                            bt_b = false;
                            t_bt = r->t - (JD_0 - ephem->jd_ref); 
                        }else{
                            bt2 = calc_bt2_particles(Ncukf, dim, CUT_M, &r->particles[Nmc + Nukf], Wcukf, mu, S);
                            if((bt2 > b2_t) || (bt2 < -b2_t)){
                                bt_b = false;
                                t_bt = r->t - (JD_0 - ephem->jd_ref); 
                            }
                        }
                    }
                    t_last = r->t; 
                }

                if(!HZ_b && !b_b && !NED_b && !bt_b) break;

                if(r->dt > max_dt) r->dt = max_dt; 
                
                reb_simulation_step(r);
            }

            if(r->t < t_end){
                #pragma omp critical
                {
                    printf("NEA: %d, Period (d): %.15f, Frob. norm: %.15e = %.15f %.15f %.15f %.15f\n", i, T, A, t_HZ, t_b, t_NED, t_bt);
                }
            }

            // Clean up memory
            free(mu); mu = NULL;
            for(int j = 0; j < dim; j++){
                free(S[j]);
            }
            free(S); S = NULL;
            free(X); X = NULL; 
            free(Zukf); free(Wukf); Zukf = NULL; Wukf = NULL;
            free(Zcukf); free(Wcukf); Zcukf = NULL; Wcukf = NULL;
            assist_free(ax); ax = NULL; 
            reb_simulation_free(r); r = NULL; 
        }
    }

    assist_ephem_free(ephem);

    double end = omp_get_wtime();
    printf("TIME: %.15f\n", end - start);

    return 0; 
}
