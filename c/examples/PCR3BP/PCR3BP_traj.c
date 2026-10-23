#include "../../src/utils.h"
#include "../../src/models.h"

int main(int argc, char* argv[]){

    

    int dim    = 4;              
    double* x  = init_v(dim); 
    double** S  = init_m(dim, dim);
    double T    = 0; 

    // directories
    const char* R_DIR = "../../../../results/PCR3BP/traj";
    const char* M_DIR = "../../measurements/PCR3BP/traj";
    char M_FILE[256];  
    snprintf(M_FILE, sizeof(M_FILE), "%s/measurement0.txt", M_DIR); 
    read_measurement(dim, M_FILE, x, S, &T); 

    char FILE_NAME[256];  // allocate enough space
    snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/PCR3BP_rv.txt", R_DIR);
    FILE* f = fopen(FILE_NAME, "w");
    fprintf(f, "t, x [LU], y [LU], vx [LU/TU], vy [LU/TU]\n"); 
    fprintf(f, "%.3f %.15e %.15e %.15e %.15e\n", 0.0, x[0], x[1], x[2], x[3]);

    double dt = 0.001;
    double t = 0.0;
    while(t < T){

        dt = fmin(dt, T - t); 
        RK4(PCR3BP, dim, x, t, dt); t += dt; 
        fprintf(f, "%.3f %.15e %.15e %.15e %.15e\n", t, x[0], x[1], x[2], x[3]);
    }

    free(x); 
    for(int i = 0; i < dim; i++){
        free(S[i]); 
    }
    free(S); 

    return 0;
}
