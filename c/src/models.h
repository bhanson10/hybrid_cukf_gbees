#ifndef MODELS_H
#define MODELS_H

#include "utils.h"

typedef void (*fn)(double, double*, double*, double*);

/*==============================================================================
                            FUNCTION DEFINITIONS
==============================================================================*/

// Dynamics models
void PCR3BP(double t, double* x, double* xp, double* coef);
void CR3BP(double t, double* x, double* xp, double* coef);

// Measurement models
void TRUTH(double t, double* x, double* y, double* coef);
void R_RR(double t, double* x, double* y, double* coef);
void R_RR_2(double t, double* x, double* y, double* coef);
void LOS(double t, double* x, double* y, double* coef);
void LOS_2(double t, double* x, double* y, double* coef);

#endif // MODELS_H
