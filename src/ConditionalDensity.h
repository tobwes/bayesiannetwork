#ifndef CONDITIONALDENSITY_H
#define CONDITIONALDENSITY_H

#include"constants.h"
#include"structures.h"
#include"includes.h"
#include"defines.h"
#include"helpers.h"

Density* initializeConditionalDensity(Density* d, Density* density, double* fixed, bool* positions, bool* parameters, bool normalize);
ConditionalDensityParameters* initializeConditionalDensityParameters(ConditionalDensityParameters* cdp, Density* density, double* fixed, bool* positions, bool *parameters, bool normalize);
double evaluateConditionalDensity(Density* this, const double* start, const double* end);
double evaluateConditionalDensityLog(Density* this, const double* start, const double* end);
void destroyConditionalDensity(Density* this);

#endif
