#ifndef LAPACK_HELPERS_H
#define LAPACK_HELPERS_H

#include"defines.h"
#include"includes.h"
#include"structures.h"

double* rep(double d, int n, double* res);
double* diag(double d, int m, int n, double* res);
void printArray(double* d, int m, int n);

#endif
