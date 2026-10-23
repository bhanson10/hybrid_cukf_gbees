#include "../../src/utils.h"
#include "../../src/models.h"

int main(int argc, char* argv[]){

    // inputs
    fn f = PCR3BP;
    int dim = 4;          
    fn h = LOS_2;
    int y_dim = 4;
    // double LU = 238529;
    // double TU = 18913;
    // double sigma[4] = {pow(1.0E-3 / LU, 2.0), pow(1.0E-3 / LU, 2.0), pow(1.0E-7 / (LU / TU), 2.0), pow(1.0E-7 / (LU / TU), 2.0)};
    double sigma[4] = {pow(1E-5, 2.0), pow(1E-5, 2.0), pow(1E-5, 2.0), pow(1E-5, 2.0)};
    const char* M_MODEL = "PCR3BP";
    int nm = 4;
    
    // parameters
    double* x   = init_v(dim); 
    double** S  = init_m(dim, dim);
    double T0    = 0; 

    // directories
    char M_DIR[256];  
    snprintf(M_DIR, sizeof(M_DIR), "../../measurements/%s", M_MODEL);
    char M_FILE[256];  
    snprintf(M_FILE, sizeof(M_FILE), "%s/measurement0.txt", M_DIR); 
    read_measurement(dim, M_FILE, x, S, &T0);
    double T = T0;  

    // Measurement parameters 
    FILE** record_files = malloc(nm * sizeof(FILE*)); 
    char FILE_NAME[256];   

    double y[y_dim]; 
    double dt = 0.001;
    double t = 0.0;
    for(int i = 1; i < nm; i++){

        while(fabs(t - T) > 1e-4){
            dt = fmin(dt, T - t); 
            RK4(f, dim, x, t, dt); 
            t+= dt; 
        }

        print_vector(dim, x); 

        // get ground station measurement
        h(t, x, y);
        print_vector(y_dim, y);         
        
        // write measurement files
        snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/measurement%d.txt", M_DIR, i);
        record_files[i] = fopen(FILE_NAME, "w");
        fprintf(record_files[i], "Mean\n"); 
        for(int j = 0; j < y_dim; j++){
            fprintf(record_files[i], "%.15e ", y[j]); 
        }
        fprintf(record_files[i], "\n\n");
        fprintf(record_files[i], "Covariance\n"); 
        for(int j = 0; j < y_dim; j++){
            for(int k = 0; k < y_dim; k++){
                if(j == k){
                    fprintf(record_files[i], "%.15e ", sigma[j]);
                }else{
                    fprintf(record_files[i], "%.15e ", 0.0);
                }
            }
            fprintf(record_files[i], "\n"); 
        }
        fprintf(record_files[i], "\n"); 
        fprintf(record_files[i], "T (TU)\n%.15f", T0); 
        fclose(record_files[i]); 

        T += T0; 
        dt = 0.001;
    }

    free(x); 
    for(int i = 0; i < dim; i++){
        free(S[i]); 
    }
    free(S); 
    free(record_files); 


    return 0;
}
