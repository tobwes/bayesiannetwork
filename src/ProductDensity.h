#ifndef PRODUCTDENSITY_H
#define PRODUCTDENSITY_H

#include"constants.h"
#include"includes.h"
#include"structures.h"
#include"defines.h"
#include"helpers.h"

Density* initializeProductDensity(Density* d, Density* d1, Density* d2);
ProductDensityParameters* initializeProductDensityParameters(ProductDensityParameters* pdp, Density* d1, Density* d2);
double evaluateProductDensity(Density* this, const double* start, const double* end);
double evaluateProductDensityLog(Density* this, const double* start, const double* end);
void destroyProductDensity(Density* this);

#endif
