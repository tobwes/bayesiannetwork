#ifndef ERRORDENSITY_H
#define ERRORDENSITY_H

#include "includes.h"
#include"defines.h"
#include"structures.h"
#include"error.h"
#include"helpers.h"

Density* initializeErrorDensity(Density* d, int type, int dimIN, int dimTARGET, int dimWeights, double beta);
ErrorDensityParameters* initializeErrorDensityParameters(ErrorDensityParameters* edp, int type, int dimIN, int dimTARGET, int dimWeights, double beta);
void destroyErrorDensity(Density* d);
double evaluateErrorDensity(Density* this, const double* first, const double* last);
double evaluateErrorDensityLog(Density* this, const double* first, const double* last);
void set_beta(Density* ed, double beta);

#endif
