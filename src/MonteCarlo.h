#ifndef MONTECARLO_H
#define MONTECARLO_H

#include"includes.h"
#include"defines.h"
#include"constants.h"
#include"structures.h"
#include"ProbabilitySpace.h"
#include"helpers.h"

MonteCarloEngine* initializeMonteCarloEngine(ProbabilitySpace* ps, Density* density, Density* sampling_prototype, int sampling_algorithm, int burnin, int gap, double rejection_sampling_factor);
void destroyMonteCarloEngine(MonteCarloEngine* mce);
const double* getRVs(ProbabilitySpace* ps, int n, bool generate_new);
void clear(MonteCarloEngine* mce);
void drawRVs(ProbabilitySpace* ps, int n); //private
void burnin(ProbabilitySpace* ps, double* start, bool new_starting_point); //private
void new_burnin(ProbabilitySpace* ps, bool new_starting_point);
void start_sampling(ProbabilitySpace* ps, int m, double* start); //private
void importance_sampling(ProbabilitySpace* ps, int m, double* start); //private
void rejection_sampling(ProbabilitySpace* ps, int m, double* start); //private
void metropolis_sampling(ProbabilitySpace* ps, int m, double* start); //private


#endif
