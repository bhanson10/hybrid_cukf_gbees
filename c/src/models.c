// models.c
// Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

#include "models.h"

void PCR3BP(double t, double* x, double* xp, double* coef){

    double mu = coef[0];
    double r1 = pow(pow(x[0] + mu, 2) + pow(x[1], 2), 1.5);
    double r2 = pow(pow(x[0] - 1 + mu, 2) + pow(x[1], 2), 1.5);

    xp[0] = x[2];
    xp[1] = x[3];
    xp[2] = 2 * x[3] + x[0] - (mu * (x[0] - 1 + mu) / r2) - ((1 - mu) * (x[0] + mu) / r1);
    xp[3] = -2 * x[2] + x[1] - (mu * x[1] / r2) - ((1 - mu) * x[1] / r1);

    return;
}

void CR3BP(double t, double* x, double* xp, double* coef){

    double mu = coef[0];
    double r1 = pow(pow(x[0] + mu, 2) + pow(x[1], 2) + pow(x[2], 2), 1.5);
    double r2 = pow(pow(x[0] - 1 + mu, 2) + pow(x[1], 2) + pow(x[2], 2), 1.5);

    xp[0] = x[3];
    xp[1] = x[4];
    xp[2] = x[5];
    xp[3] = 2 * x[4] + x[0] - (mu * (x[0] - 1 + mu) / r2) - ((1 - mu) * (x[0] + mu) / r1);
    xp[4] = -2 * x[3] + x[1] - (mu * x[1] / r2) - ((1 - mu) * x[1] / r1);
    xp[5] = -(mu * x[2] / r2) - ((1 - mu) * x[2] / r1);

    return;
}

void TRUTH(double t, double* x, double* y, double* coef){

    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
    y[3] = x[3];

    return;
}

void R_RR(double t, double* x, double* y, double* coef){

    double x_a[6] = {0.83691513, 0.0, 0.0, 0.0, 0.0, 0.0}; 
    double* rho = init_v(3); 
    double* rho_dot = init_v(3); 
    for(int i = 0; i < 3; i++){
        rho[i] = x[i] - x_a[i]; 
        rho_dot[i] = x[i + 3] - x_a[i]; 
    }

    y[0] = sqrt(rho[0] * rho[0] + rho[1] * rho[1] + rho[2] * rho[2]); // range
    y[1] = dot(rho, rho_dot, 3) / y[0]; // range-rate

    free(rho); 
    free(rho_dot); 

    return;
}

void R_RR_2(double t, double* x, double* y, double* coef){

    double L1 = 0.99601828;
    double dx1 = x[0] - L1;
    double L2 = 1.00399194;
    double dx2 = x[0] - L2;
    double dy = x[1];

    y[0] = sqrt(dx1 * dx1 + dy * dy);
    y[1] = sqrt(dx2 * dx2 + dy * dy);
    y[2] = (dx1 * x[2] + dy * x[3]) / y[0];
    y[3] = (dx2 * x[2] + dy * x[3]) / y[1];

    return; 
}

void LOS(double t, double* x, double* y, double* coef){

    double L1 = 0.99601828;
    double dx = x[0] - L1;
    double dy = x[1];
    double rho = sqrt(dx * dx + dy * dy);

    y[0] = dx / rho; 
    y[1] = dy / rho; 

    return; 
}

void LOS_2(double t, double* x, double* y, double* coef){

    double L1 = 0.99601828;
    double dx1 = x[0] - L1;
    double L2 = 1.00399194;
    double dx2 = x[0] - L2;
    double dy = x[1];
    double rho1 = sqrt(dx1 * dx1 + dy * dy);
    double rho2 = sqrt(dx2 * dx2 + dy * dy);

    y[0] = dx1 / rho1; 
    y[1] = dy / rho1; 
    y[2] = dx2 / rho2; 
    y[3] = dy / rho2; 

    return; 
}