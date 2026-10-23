// pme.c
// Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

#include "pme.h"

static inline double ipow(double x, int e){
    double r = 1.0;
    for (int i = 0; i < e; i++) r *= x;
    return r;
}

void dual_objective(const double* l, int p, int N, double* Phi, double* Mi, double w, double* f, double* g, double* H, double* s, double* e_scaled, double* Phi_scaled, bool grad, bool hess){

    // compute s = Phi * l and find smax simultaneously
    double smax = -HUGE_VAL;
    for(int i = 0; i < N; i++){
        const double* Phi_i = Phi + i * p;
        double si = 0.0;
        for(int j = 0; j < p; j++)
            si += Phi_i[j] * l[j];
        s[i] = si;
        if(si > smax)
            smax = si;
    }

    // compute partition function
    double scale = w * exp(smax);
    double Z = 0.0;
    for(int i = 0; i < N; i++){
        double ei = exp(s[i] - smax);
        e_scaled[i] = scale * ei;
        Z += ei;
    }

    // compute objective
    double lM = 0.0;
    for(int i = 0; i < p; i++)
        lM += l[i] * Mi[i];
    *f = scale * Z - lM;

    // compute gradient 
    if(grad){

        for(int i = 0; i < p; i++) g[i] = -Mi[i];
        cblas_dgemv(CblasRowMajor, CblasTrans, N, p, 1.0, Phi, p, e_scaled, 1, 1.0, g, 1);
    }

    // compute Hessian 
    if(hess){
        for(int i = 0; i < N; i++){

            const double weight = e_scaled[i];
            const double* Phi_i = Phi + i * p;
            double* Phi_scaled_i = Phi_scaled + i * p;

            for(int j = 0; j < p; j++)
                Phi_scaled_i[j] = weight * Phi_i[j];
        }

        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, p, p, N, 1.0, Phi, p, Phi_scaled, p, 0.0, H, p);
    }

    return; 
}

double comp_count(int r, int k){
    return nchoosek(r + k - 1, k);
}

void gen_E(int d, int M, double idx, int* vec){

    double remaining_idx = idx;
    for(int k = 0; k <= M; k++){

        double block_size = comp_count(d, k);
        if(remaining_idx > block_size){
            remaining_idx -= block_size;
            continue;
        }

        int rem_dim = d;
        int rem_sum = k;
        for(int i = 0; i < d - 1; i++){

            for(int val = rem_sum; val >= 0; val--){
                double count = comp_count(rem_dim - 1, rem_sum - val);
                if(remaining_idx > count){
                    remaining_idx -= count;
                }else{
                    vec[i] = val;
                    rem_sum -= val;
                    rem_dim--;
                    break;
                }
            }
        }

        vec[d - 1] = rem_sum;
        return;
    }
}

void cut2mom(double** Z, double* W, int N, int M, int dim, int p, double* Mi, int** E){

    for(int i = 0; i < p; i++){
        Mi[i] = 0.0; 

        gen_E(dim, M, i + 1, E[i]); 

        for(int j = 0; j < N; j++){
            double f = 1.0; 
            for(int k = 0; k < dim; k++){
                int e = E[i][k];
                if (e)
                    f *= ipow(Z[j][k], e);
            }
            Mi[i] += W[j] * f; 
        }
    }
    
    return; 
}

double* solve_pme_eqns(double* Mi, int** E, double* X, double vol, int N, int dim, int p){

    // ============================================================
    // Gaussian quadrature
    // ============================================================
    double w = vol / N;

    // ============================================================
    // Build basis matrix
    // ============================================================
    double* Phi = init_v(N * p);
    for(int k = 0; k < N; k++){
        const double* Xk = X + k * dim;
        double* Phi_k = Phi + k * p;
        for(int i = 0; i < p; i++){
            double val = 1.0;
            for(int j = 0; j < dim; j++){
                int e = E[i][j];
                if(e != 0)
                    val *= ipow(Xk[j], e);
            }
            Phi_k[i] = val;
        }
    }

    // ============================================================
    // Newton settings
    // ============================================================
    const int    max_iter   = 100;           // Gradient convergence tolerance
    const double grad_tol   = 1e-8;
    const double step_tol   = 1e-10;         // Step convergence tolerance
    const double f_tol      = 1e-12;         // Objective convergence tolerance
    const double stagnation_grad_tol = 1e-6; // Gradient tolerance used when objective has stagnated
    const double beta       = 0.5;           // Backtracking parameters
    const double c1         = 1e-4;
    const double min_step   = 1e-12;

    // ============================================================
    // Allocate workspace
    // ============================================================
    double* l       = init_v(p);
    double* g       = init_v(p);
    double* H       = init_v(p * p);
    double* R       = init_v(p * p);
    double* delta   = init_v(p);
    double* H_reg   = init_v(p * p);
    double* l_trial = init_v(p);
    for(int i = 0; i < p; i++) l[i] = 0.0;

    // ============================================================
    // Workspace for dual_objective
    // ============================================================
    double* s           = init_v(N);
    double* e_scaled    = init_v(N);
    double* Phi_scaled  = init_v(N * p);

    // ============================================================
    // Display
    // ============================================================
    bool display_flag = true;
    if(display_flag){
        printf("\n");
        printf("PME Newton Solver\n");
        printf("-----------------\n");
        printf("%5s %16s %16s %16s %16s %10s\n",
               "Iter",
               "Objective",
               "|dF|",
               "||Gradient||",
               "||Step||",
               "Alpha");
    }


    // ============================================================
    // Newton iteration
    // ============================================================
    double f  = 0.0;
    double fl = 0.0;
    for(int i = 1; i <= max_iter; i++){

        // --------------------------------------------------------
        // Evaluate objective, gradient, and Hessian
        // --------------------------------------------------------
        dual_objective(l, p, N, Phi, Mi, w, &f, g, H, s, e_scaled, Phi_scaled, true, true);
        double grad_norm = norm_inf(g, p);
        double f_change = 0.0;
        if(i > 1) f_change = fabs(f - fl);

        // --------------------------------------------------------
        // Convergence: gradient
        // --------------------------------------------------------
        if(grad_norm < grad_tol){
            if(display_flag){
                printf("%5d %16.8e %16.8e %16.8e %16s %10s\n",
                       i - 1,
                       f,
                       f_change,
                       grad_norm,
                       "---",
                       "---");
                printf("Converged: gradient tolerance.\n");
            }

            free(Phi);
            free(g);
            free(H);
            free(R);
            free(delta);
            free(H_reg);
            free(l_trial);
            free(s);
            free(e_scaled);
            free(Phi_scaled);

            return l;
        }


        // --------------------------------------------------------
        // Convergence: objective stagnation
        //
        // Only use this if the gradient is already reasonably
        // small. This prevents declaring convergence simply because
        // the objective is flat far away from the solution.
        // --------------------------------------------------------
        if(i > 1){
            double f_scale = 1.0 + fabs(f);
            if(f_change < f_tol * f_scale && grad_norm < stagnation_grad_tol){
                if(display_flag){
                    printf("%5d %16.8e %16.8e %16.8e %16s %10s\n",
                           i - 1,
                           f,
                           f_change,
                           grad_norm,
                           "---",
                           "---");

                    printf("Converged: objective stagnation.\n");
                }

                free(Phi);
                free(g);
                free(H);
                free(R);
                free(delta);
                free(H_reg);
                free(l_trial);
                free(s);
                free(e_scaled);
                free(Phi_scaled);

                return l;
            }
        }

        // --------------------------------------------------------
        // Symmetrize Hessian
        // --------------------------------------------------------
        for(int j = 0; j < p; j++){
            for(int k = 0; k <= j; k++){
                double sym_val = 0.5 * (H[j * p + k] + H[k * p + j]);
                H[j * p + k] = sym_val;
                H[k * p + j] = sym_val;
            }
        }

        // --------------------------------------------------------
        // Cholesky decomposition
        // --------------------------------------------------------
        int flag = chol2(p, H, R);
        if(flag == 0){

            // ----------------------------------------------------
            // delta = -H^{-1} g
            //
            // H = R R^T
            // ----------------------------------------------------

            for(int j = 0; j < p; j++)
                delta[j] = -g[j];

            cblas_dtrsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, p, R, p, delta, 1);

            cblas_dtrsv(CblasRowMajor, CblasLower, CblasTrans, CblasNonUnit, p, R, p, delta, 1);

        }else{

            // ----------------------------------------------------
            // Hessian regularization
            // ----------------------------------------------------

            double regularization =
                1e-10 * fmax(1.0, norm_frob(H, p * p));

            for(int j = 0; j < p; j++){

                for(int k = 0; k < p; k++){

                    H_reg[j * p + k] = H[j * p + k];

                    if(j == k)
                        H_reg[j * p + k] += regularization;
                }
            }


            // ----------------------------------------------------
            // Cholesky of regularized Hessian
            // ----------------------------------------------------

            flag = chol2(p, H_reg, R);


            if(flag == 0){

                for(int j = 0; j < p; j++)
                    delta[j] = -g[j];

                cblas_dtrsv(
                    CblasRowMajor,
                    CblasLower,
                    CblasNoTrans,
                    CblasNonUnit,
                    p,
                    R,
                    p,
                    delta,
                    1
                );

                cblas_dtrsv(
                    CblasRowMajor,
                    CblasLower,
                    CblasTrans,
                    CblasNonUnit,
                    p,
                    R,
                    p,
                    delta,
                    1
                );

            }else{

                fprintf(
                    stderr,
                    "solve_pme_eqns: Newton Hessian is not positive definite. "
                    "Pseudoinverse required.\n"
                );

                free(Phi);
                free(g);
                free(H);
                free(R);
                free(delta);
                free(H_reg);
                free(l_trial);
                free(s);
                free(e_scaled);
                free(Phi_scaled);

                return NULL;
            }
        }


        // --------------------------------------------------------
        // Check Newton step
        // --------------------------------------------------------

        double step_norm = norm_inf(delta, p);

        double l_norm = norm_inf(l, p);

        if(step_norm < step_tol * (1.0 + l_norm)){

            // Take the Newton step
            for(int j = 0; j < p; j++)
                l[j] += delta[j];


            if(display_flag){

                printf("%5d %16.8e %16.8e %16.8e %16.8e %10.4f\n",
                       i - 1,
                       f,
                       f_change,
                       grad_norm,
                       step_norm,
                       1.0);

                printf("Converged: step tolerance.\n");
            }


            free(Phi);
            free(g);
            free(H);
            free(R);
            free(delta);
            free(H_reg);
            free(l_trial);
            free(s);
            free(e_scaled);
            free(Phi_scaled);

            return l;
        }


        // --------------------------------------------------------
        // Check whether Newton direction is a descent direction
        // --------------------------------------------------------

        double directional_derivative = 0.0;

        for(int j = 0; j < p; j++)
            directional_derivative += g[j] * delta[j];


        // If numerical issues have caused the Newton direction to
        // stop being a descent direction, do not blindly perform
        // an Armijo line search.
        if(directional_derivative >= 0.0){

            fprintf(
                stderr,
                "solve_pme_eqns: Newton direction is not a descent "
                "direction (g^T delta = %.6e).\n",
                directional_derivative
            );

            // If we are already very close to a stationary point,
            // accept this as convergence.
            if(grad_norm < stagnation_grad_tol){

                if(display_flag)
                    printf("Converged: near-stationary point.\n");

                free(Phi);
                free(g);
                free(H);
                free(R);
                free(delta);
                free(H_reg);
                free(l_trial);
                free(s);
                free(e_scaled);
                free(Phi_scaled);

                return l;
            }

            free(Phi);
            free(g);
            free(H);
            free(R);
            free(delta);
            free(H_reg);
            free(l_trial);
            free(s);
            free(e_scaled);
            free(Phi_scaled);

            return NULL;
        }


        // --------------------------------------------------------
        // Backtracking line search
        // --------------------------------------------------------

        double alpha = 1.0;

        bool accepted = false;

        while(alpha > min_step){

            // ----------------------------------------------------
            // Trial point
            // ----------------------------------------------------

            for(int j = 0; j < p; j++)
                l_trial[j] = l[j] + alpha * delta[j];


            // ----------------------------------------------------
            // Evaluate trial objective
            // ----------------------------------------------------

            double f_trial;

            dual_objective(
                l_trial,
                p,
                N,
                Phi,
                Mi,
                w,
                &f_trial,
                NULL,
                NULL,
                s,
                e_scaled,
                Phi_scaled,
                false,
                false
            );


            // ----------------------------------------------------
            // Armijo condition
            // ----------------------------------------------------

            double armijo_term =
                c1 * alpha * directional_derivative;

            double armijo_rhs = f + armijo_term;


            // Standard Armijo condition
            if(f_trial <= armijo_rhs){

                accepted = true;
                break;
            }


            // ----------------------------------------------------
            // Floating-point safeguard
            //
            // If the required decrease is smaller than the
            // resolution of double precision relative to f, then
            // further backtracking is meaningless.
            // ----------------------------------------------------

            double machine_resolution =
                10.0 * DBL_EPSILON * (1.0 + fabs(f));

            if(fabs(armijo_term) <= machine_resolution){

                accepted = true;
                break;
            }


            // ----------------------------------------------------
            // Backtrack
            // ----------------------------------------------------

            alpha *= beta;
        }


        // --------------------------------------------------------
        // Line-search failure
        // --------------------------------------------------------

        if(!accepted){

            // ----------------------------------------------------
            // If the objective is already numerically stagnant and
            // the gradient is reasonably small, consider this
            // convergence rather than failure.
            // ----------------------------------------------------

            double f_scale = 1.0 + fabs(f);

            if(grad_norm < stagnation_grad_tol &&
               f_change < f_tol * f_scale){

                if(display_flag){

                    printf("%5d %16.8e %16.8e %16.8e %16.8e %10.4e\n",
                           i - 1,
                           f,
                           f_change,
                           grad_norm,
                           alpha * step_norm,
                           alpha);

                    printf(
                        "Converged: line search reached "
                        "floating-point precision.\n"
                    );
                }

                free(Phi);
                free(g);
                free(H);
                free(R);
                free(delta);
                free(H_reg);
                free(l_trial);
                free(s);
                free(e_scaled);
                free(Phi_scaled);

                return l;
            }


            // ----------------------------------------------------
            // Genuine line-search failure
            // ----------------------------------------------------

            fprintf(
                stderr,
                "solve_pme_eqns: Newton line search failed.\n"
            );

            if(display_flag){

                printf("%5d %16.8e %16.8e %16.8e %16.8e %10.4e\n",
                       i - 1,
                       f,
                       f_change,
                       grad_norm,
                       alpha * step_norm,
                       alpha);
            }


            free(Phi);
            free(g);
            free(H);
            free(R);
            free(delta);
            free(H_reg);
            free(l_trial);
            free(s);
            free(e_scaled);
            free(Phi_scaled);

            return l;
        }


        // --------------------------------------------------------
        // Newton update
        // --------------------------------------------------------

        for(int j = 0; j < p; j++)
            l[j] += alpha * delta[j];


        // --------------------------------------------------------
        // Display
        // --------------------------------------------------------

        if(display_flag){

            printf("%5d %16.8e %16.8e %16.8e %16.8e %10.4f\n",
                   i - 1,
                   f,
                   f_change,
                   grad_norm,
                   alpha * step_norm,
                   alpha);
        }


        // --------------------------------------------------------
        // Store objective for next iteration
        // --------------------------------------------------------

        fl = f;
    }


    // ============================================================
    // Maximum iterations
    // ============================================================

    fprintf(
        stderr,
        "solve_pme_eqns: Newton solver reached maximum iterations.\n"
    );


    // ============================================================
    // Cleanup
    // ============================================================

    free(Phi);
    free(g);
    free(H);
    free(R);
    free(delta);
    free(H_reg);
    free(l_trial);
    free(s);
    free(e_scaled);
    free(Phi_scaled);


    return l;
}

int whiten_cut(double** Z, double* W, int N, int M, int dim, double** Zp, double* mu, double* L){

    // compute mean
    double* mut = reconstruct_cut(Z, W, dim, 1, M);
    for(int i = 0; i < dim; i++) mu[i] = mut[i]; 

    // compute covariance
    double* S = reconstruct_cut(Z, W, dim, 2, M);

    // cholesky decomposition
    int flag = chol2(dim, S, L);
    if(flag != 0) return -2;

    // center sigma points
    double* Zpt = init_v(N * dim); 
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < dim; ++j){
            Zpt[i * dim + j] = Z[i][j] - mu[j];
        }
    }

    // solve for whitened cut points
    for (int i = 0; i < N; ++i){
        cblas_dtrsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, dim, L, dim, &Zpt[i * dim], 1);
    }
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < dim; ++j){
            Zp[i][j] = Zpt[i * dim + j];
        }
    }

    free(S);

    return 0;
}

double* phi(double* X, int N, double* l, int** E, int p, int dim){
    double* f = init_v(N);
    if (!f) return NULL;

    const double max_mahalanobis_dist = 4.0;
    const double max_dist_sq = max_mahalanobis_dist * max_mahalanobis_dist;

    double max_fk = -DBL_MAX;

    // Pass 1: Compute log-densities (fk) and locate the peak exponent max_fk
    for (int k = 0; k < N; k++) {
        double* Xk = X + k * dim;

        // Distance check
        double dist_sq = 0.0;
        for (int j = 0; j < dim; j++) {
            dist_sq += Xk[j] * Xk[j];
        }

        if (dist_sq > max_dist_sq) {
            f[k] = -DBL_MAX; // Outside support
            continue;
        }

        // Polynomial log-density evaluation
        double fk = 0.0;
        for (int i = 0; i < p; i++) {
            double gi = 1.0;
            for (int j = 0; j < dim; j++) {
                int e = E[i][j];
                if (e == 1) {
                    gi *= Xk[j];
                } else if (e > 1) {
                    gi *= ipow(Xk[j], e);
                }
            }
            fk += l[i] * gi;
        }

        f[k] = fk;
        if (fk > max_fk) {
            max_fk = fk;
        }
    }

    // Pass 2: Exponentiate using max-shift offset to guarantee values in [0, 1]
    for (int k = 0; k < N; k++) {
        if (f[k] == -DBL_MAX) {
            f[k] = 0.0;
            continue;
        }

        double shifted_fk = f[k] - max_fk; // Guarantee shifted_fk <= 0

        if (shifted_fk < -700.0) {
            f[k] = 0.0; // Underflow in the deep tails
        } else {
            f[k] = exp(shifted_fk); // Core points will now be near 1.0
        }
    }

    return f;
}

void pme_cut2pdf(double** Z, double* W, int N, int CUT_M, int dim, int PME_M, int n_quad, int n_eval, double** X, double** P){

    // size
    int num_pts_quad = 1;
    for (int i = 0; i < dim; i++) num_pts_quad *= n_quad;
    int num_pts_eval = 1;
    for (int i = 0; i < dim; i++) num_pts_eval *= n_eval;

    // numerical preconditioning
    double** Zp = init_m(N, dim); 
    double* mu = init_v(dim); 
    double* L = init_v(dim * dim); 
    whiten_cut(Z, W, N, CUT_M, dim, Zp, mu, L);

    // converting CUT to scalar moment values
    int p = nchoosek(dim + PME_M, dim);
    double* Mi = init_v(p); 
    int** E = init_m_i(p, dim); 
    cut2mom(Zp, W, N, PME_M, dim, p, Mi, E);

    printf("\nMi=\n");
    print_vector_T(p, Mi); 
    printf("E=\n");
    print_matrix_i(p, dim, E); 

    // -------------------------------------------------------
    // solving for lagrange multipliers
    // -------------------------------------------------------

    // constructing support for Gaussian quadrature integration
    double** Omega = init_m(dim, 2); 
    for(int i = 0; i < dim; i++){
        double Zmin = Zp[0][i]; double Zmax = Zmin; 
        for(int j = 1; j < N; j++){
            Zmin = fmin(Zmin, Zp[j][i]); 
            Zmax = fmax(Zmax, Zp[j][i]); 
        }
        Omega[i][0] = Zmin - 0.1 * (Zmax - Zmin); 
        Omega[i][1] = Zmax + 0.1 * (Zmax - Zmin); 
    }
    double vol = 1.0;
    for(int i = 0; i < dim; i++){
        vol *= Omega[i][1] - Omega[i][0];
    }

    // constructing Gaussian quadrature points
    double* Xp = init_v(fmax(num_pts_quad, num_pts_eval) * dim); 
    for(int i = 0; i < num_pts_quad; i++){
        int idx = i;

        for(int j = 0; j < dim; j++){
            int digit = idx % n_quad;
            idx /= n_quad;
            double h = (n_quad > 1) ? (Omega[j][1] - Omega[j][0]) / (n_quad - 1) : 0.0;
            Xp[i * dim + j] = Omega[j][0] + digit * h;
        }
    }

    double* l = solve_pme_eqns(Mi, E, Xp, vol, num_pts_quad, dim, p); 
    printf("\n"); 
    printf("l=\n");
    print_vector_T(p, l); 

    // construct evaluation points in original space
    for(int i = 0; i < dim; i++){
        double Zmin = Z[0][i]; double Zmax = Zmin; 
        for(int j = 1; j < N; j++){
            Zmin = fmin(Zmin, Z[j][i]); 
            Zmax = fmax(Zmax, Z[j][i]); 
        }
        Omega[i][0] = Zmin - 0.2 * (Zmax - Zmin); 
        Omega[i][1] = Zmax + 0.2 * (Zmax - Zmin); 
    }
    
    *X = init_v(num_pts_eval * dim); 
    for(int i = 0; i < num_pts_eval; i++){
        int idx = i;

        for(int j = 0; j < dim; j++){
            int digit = idx % n_eval;
            idx /= n_eval;
            double h = (n_eval > 1) ? (Omega[j][1] - Omega[j][0]) / (n_eval - 1) : 0.0;
            (*X)[i * dim + j] = Omega[j][0] + digit * h;
        }
    }

    // whiten evaluation points
    for (int i = 0; i < num_pts_eval; ++i) {
        for (int j = 0; j < dim; ++j) {
            Xp[i * dim + j] = (*X)[i * dim + j] - mu[j];
        }
        cblas_dtrsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, dim, L, dim, &Xp[i * dim], 1);
    }

    // calculate PDF over evaluation points
    *P = phi(Xp, num_pts_eval, l, E, p, dim);
    double sum = 0.0;
    for(int i = 0; i < num_pts_eval; i++) sum += (*P)[i];
    for(int i = 0; i < num_pts_eval; i++) (*P)[i] /= sum;
    
    // cleanup
    for(int i = 0; i < dim; i++) free(Omega[i]); 
    free(Omega); 
    free(Xp);
    free(Zp);
    free(mu);
    free(L);
    free(Mi);

    return; 
}
