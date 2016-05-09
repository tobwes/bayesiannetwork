#ifndef GNEURAL_NETWORK_H
#define GNEURAL_NETWORK_H

#include"ProbabilitySpace.h"
#include"GaussianFamily.h"
#include"ErrorDensity.h"
#include"ConditionalDensity.h"
#include"ProductDensity.h"
#include"error.h"
#include"lapack_helpers.h"


void evaluateNetworkBayesian();
void outputRV(const double* start, const double* end, double* y);

#endif
