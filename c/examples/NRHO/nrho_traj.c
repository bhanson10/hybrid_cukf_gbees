#include "../../src/utils.h"
#include "../../src/models.h"

int main(int argc, char* argv[]){

    // inputs
    fn f = CR3BP;
    int dim = 6;         
    const char* orbit = "NRHO"; 
    const char* len = "30d"; 
    double coef_f[1] = {0.012150585609624}; // Earth-Moon system

    // parameters       
    double* x  = init_v(dim); 
    double** S = init_m(dim, dim);
    double T   = 0; 
    char M_FILE[256];  
    snprintf(M_FILE, sizeof(M_FILE), "../../measurements/%s/traj/measurement0_%s.txt", orbit, len); 
    read_measurement(dim, M_FILE, x, S, &T); 

    // directories
    char R_DIR[256];  
    snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s", orbit);
    mkdir(R_DIR, 0777);
    snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s/traj", orbit);
    mkdir(R_DIR, 0777);
    snprintf(R_DIR, sizeof(R_DIR), "../../../../results/%s/traj/%s", orbit, len);
    mkdir(R_DIR, 0777);

    char FILE_NAME[256];  // allocate enough space
    snprintf(FILE_NAME, sizeof(FILE_NAME), "%s/%s_rv.txt", R_DIR, orbit);
    FILE* file = fopen(FILE_NAME, "w");
    fprintf(file, "t, x [LU], y [LU], z [LU], vx [LU/TU], vy [LU/TU], vz [LU/TU]\n"); 
    fprintf(file, "%.3f %.15e %.15e %.15e %.15e %.15e %.15e\n", 0.0, x[0], x[1], x[2], x[3], x[4], x[5]);

    double dt = 0.001;
    double t = 0.0;
    while(t < T){
        dt = fmin(dt, T - t); 
        RK4(f, dim, x, coef_f, t, dt); t += dt; 
        fprintf(file, "%.3f %.15e %.15e %.15e %.15e %.15e %.15e\n", t, x[0], x[1], x[2], x[3], x[4], x[5]);
    }

    free(x); 
    for(int i = 0; i < dim; i++){
        free(S[i]); 
    }
    free(S); 

    return 0;
}
