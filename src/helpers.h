#ifndef HELPERS_H
#define HELPERS_H

#include"includes.h"

void vectorToString(const double* vec, int dim, char* str);

void mean(const double* rvs, int dim, int nrvs, double* m);
void variance(const double* rvs, int dim, int nrvs, const double *m, double *var);
void quantile(const double* rvs, int nrvs, double *q, int n_quantiles, double* res);

int double_leq(const void* p1, const void* p2);

#endif
