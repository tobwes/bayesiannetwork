#include"includes.h"
#include"defines.h"
#include"ProbabilitySpace.h"
#include"Functions.h"

#ifndef EXPONENTIALFAMILY_H
#define EXPONENTIALFAMILY_H

/** _eta will be _copied_ into this struct **/
/** theta will _not_ be copied into this struct, it should, therefore, be allocated using malloc and will be free'd by destroyExponentialParameters unless otherwise specified (by setting free_theta = false) **/
ExponentialParameters* initializeExponentialParameters(ExponentialParameters* ep,
						       int dim,
						       int dimTheta,
						       int dimScalarProduct,
						       double* theta,
						       double _g,
						       double _A,
						       double *_eta,
						       void (*g)(const double* first, const double* last, double* res),
						       void (*A)(const double* first, const double* last, double* res),
						       void (*eta)(const double* first, const double* last, double* res),
						       void (*h)(const double* first, const double* last, double* res),
						       void (*T)(const double* first, const double* last, double* res));
void destroyExponentialParameters(ExponentialParameters* ep);

/** _eta will be _copied_ into this struct **/
/** theta will _not_ be copied into this struct, it should, therefore, be allocated using malloc and will be free'd by destroyExponentialFamily unless otherwise specified (by setting free_theta = false) **/
Density* initializeExponentialFamily(Density* ef,
				     int dim,
				     int dimTheta,
				     int dimScalarProduct,
				     double* theta,
				     double _g,
				     double _A,
				     double *_eta,
				     void (*g)(const double* first, const double* last, double* res),
				     void (*A)(const double* first, const double* last, double* res),
				     void (*eta)(const double* first, const double* last, double* res),
				     void (*h)(const double* first, const double* last, double* res),
				     void (*T)(const double* first, const double* last, double* res));
void destroyExponentialFamily(Density* ef);

Density* initializeExponential(Density* ed, double lambda);
void destroyExponential(Density* e);

double evaluateExponentialFamily(Density* this, const double* first, const double* last);
double evaluateExponentialFamilyLog(Density* this, const double* first, const double* last);

#endif
