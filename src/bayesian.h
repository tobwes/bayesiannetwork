#ifndef BAYESIAN_H
#define BAYESIAN_H

#include"includes.h"
#include"defines.h"
#include"structures.h"
#include"ProbabilitySpace.h"
#include"GaussianFamily.h"
#include"ErrorDensity.h"
#include"ConditionalDensity.h"
#include"ProductDensity.h"
#include"error.h"
#include"lapack_helpers.h"
#include"helpers.h"

extern ProbabilitySpace PROBABILITYSPACE;
extern int NDATA;
extern int NNUM;
extern int ERROR_TYPE;
extern int MAX_OPTIMIZATION_RUNS;
extern double ERROR_BOUND;
extern double PROTOTYPE_VARIANCE;

void setup_probabilityspace(int verbosity);
double calculate_training_error();
void errorFunction(double* res);

#endif
