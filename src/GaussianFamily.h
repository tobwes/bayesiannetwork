#include"includes.h"
#include"defines.h"
#include"constants.h"
#include"structures.h"
#include"ProbabilitySpace.h"

#include<lapacke.h>

#ifndef GAUSSIANFAMILY_H
#define GAUSSIANFAMILY_H

/** sigma must be allocated with malloc and will be free'd by calling destroyGaussianVarianceMatrix **/
GaussianVarianceMatrix* initializeGaussianVarianceMatrix(GaussianVarianceMatrix* gvm, int dim, double *sigma, int factorized);
void destroyGaussianVarianceMatrix(GaussianVarianceMatrix* gvm);

/** sigma must be allocated with malloc and will be free'd by calling destroyGaussianParameters **/
/** mu should be allocated with malloc and will be free'd by calling destroyGaussianParameters, as long as this is not changed manually (by setting free_mu = false) **/
GaussianParameters* initializeGaussianParameters(GaussianParameters* gp, int dim, const double* mu, double* sigma);
void destroyGaussianParameters(GaussianParameters* gp);

/** sigma must be allocated with malloc and will be free'd by calling destroyGaussian **/
/** mu should be allocated with malloc and will be free'd by calling destroyGaussian, as long as this is not changed manually (by setting free_mu = false) **/
Density* initializeGaussian(Density* d, int dim, const double* mu, double* sigma);
void destroyGaussian(Density* d);

void solveGaussianLinearEquation(GaussianVarianceMatrix* gvm, double* rhs, int NRHS, double* x);

double evaluateGaussian(Density* this, const double* first, const double* last);
double evaluateGaussianLog(Density* this, const double* first, const double* last);

void gaussianTransition(Density* d, const double* first, const double* last); // takes ownership of first - last...
double* gaussianRnd(Density* d, double* first);

#endif
