#include "../../src/utils.h"
#include "../../src/models.h"
#include "../../src/cukf.h"
#include "../../src/pme.h"
#include <gsl/gsl_cdf.h>

int main(int argc, char* argv[]){
    
    int nm = 3; 
    int dim = 4; 
    double** X = NULL;
    double* P = NULL;
    int N = 0;             
    
    // directories
    const char* M_DIR = "../../../../results/PCR3BP/cukf_gbees/posteriors/";
    char M_FILE[256];  
    snprintf(M_FILE, sizeof(M_FILE), "%s/P%d_pdf_0.bin", M_DIR, nm); 
    parse_gbees(M_FILE, dim, &P, &X, &N);

    double* mu = init_v(dim);
    double** S = init_m(dim, dim);
    weighted_mean_cov(X, P, N, dim, mu, S);

    print_vector(dim, mu);
    print_matrix(dim, dim, S);

    // free memory
    free(mu); 
    for(int i = 0; i < dim; i++){
        free(S[i]); 
    }
    free(S); 
    for(int i = 0; i < N; i++){
        free(X[i]); 
    }
    free(X); free(P); 
 
    return 0; 
}
