#ifndef PROBABILITYSPACE_H
#define PROBABILITYSPACE_H

#include"includes.h"
#include"defines.h"
#include"structures.h"
#include"MonteCarlo.h"
#include"network_helpers.h"

ProbabilitySpace* initializeProbabilitySpace(ProbabilitySpace* ps, int dim, Density* d);
void destroyProbabilitySpace(ProbabilitySpace* ps);
void registerRVs(ProbabilitySpace* ps, double* begin, int n);
void registerRV(ProbabilitySpace* ps, double* rv);
void createRVs(ProbabilitySpace* ps);
const double* getIID(ProbabilitySpace* ps, int n, bool generate_new, const double* res); // a copy of the last state is also saved internally in the registered RVs
const double* getRand(ProbabilitySpace* ps, bool generate_new);
void expectedValue(ProbabilitySpace* ps, void (*f)(double* res), int dimf, int nrvs, double* res);
void expectedValueFunctional(ProbabilitySpace* ps, void (*f)(const double* first, const double* last, double* res), int dimf, int nrvs, bool newRVs, double* res);
void setupMonteCarloEngine(ProbabilitySpace* ps, Density* sampling_prototype, int sampling_algorithm, int burnin, int gap, ...);
void setState(ProbabilitySpace* ps, double* start);



#endif
