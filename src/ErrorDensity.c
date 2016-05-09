#include"ErrorDensity.h"

ErrorDensityParameters* initializeErrorDensityParameters(ErrorDensityParameters* edp, int type, int dimIN, int dimTARGET, int dimWEIGHTS, double beta){
  if(edp==0) edp=malloc(sizeof(ErrorDensityParameters));

  edp->type = type;
  edp->beta = beta;
  edp->dimIN = dimIN;
  edp->dimTARGET = dimTARGET;
  edp->dimWEIGHTS = dimWEIGHTS;
  
  return edp;
}

Density* initializeErrorDensity(Density* d, int type, int dimIN, int dimTARGET,  int dimWEIGHTS, double beta){
  if(d==0) d=malloc(sizeof(Density));
  
  d->parameters = initializeErrorDensityParameters(0,type, dimIN, dimTARGET, dimWEIGHTS, beta);

  d->dim = dimIN + dimTARGET + dimWEIGHTS;
  d->normalizing_factor = 1; // unknown
  d->log_normalizing_factor = 0; // unknown
  
  d->evaluate = &evaluateErrorDensity;
  d->evaluateLog = &evaluateErrorDensityLog;
  d->destroy = &destroyErrorDensity;
  d->rnd = 0;
  d->setTransitionPoint = 0;
  
  return d;
}

void destroyErrorDensity(Density* d){
  ErrorDensityParameters* edp = (ErrorDensityParameters*)d->parameters;

  free(edp);
}

double evaluateErrorDensity(Density* this, const double* first, const double* last){
  ErrorDensityParameters* edp = (ErrorDensityParameters*)this->parameters;

  const double *in=first, *target=first+edp->dimIN, *weights=first+edp->dimIN+edp->dimTARGET;
  double err = specificError(edp->type,in,target,weights,false);

  return this->normalizing_factor*exp(-edp->beta*err);
}

double evaluateErrorDensityLog(Density* this, const double* first, const double* last){
  ErrorDensityParameters* edp = (ErrorDensityParameters*)this->parameters;

  const double *in=first, *target=first+edp->dimIN, *weights=target+edp->dimTARGET;
    
  double err = specificError(edp->type,in,target,weights,false);
  
  return this->log_normalizing_factor - edp->beta*err;
}

void set_beta(Density* ed, double beta){
  ErrorDensityParameters* edp = (ErrorDensityParameters*) ed->parameters;

  edp->beta = beta;
}
