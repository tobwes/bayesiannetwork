#include"ConditionalDensity.h"

// fixed needs to be a vector of length density->dim, not of the number of fixed positions!!!
Density* initializeConditionalDensity(Density* d, Density* density, double* fixed, bool* positions, bool *parameters, bool normalize){
  if(d==0) d=malloc(sizeof(Density));

  d->parameters = initializeConditionalDensityParameters(0,density,fixed,positions,parameters,normalize);
  
  d->dim = density->dim;
  for(int i=0; i<density->dim; i++){
    if(positions[i]) d->dim--;
  }
  
  d->normalizing_factor = 1;
  d->log_normalizing_factor = 0;
  
  d->evaluate = &evaluateConditionalDensity;
  if(density->evaluateLog){ d->evaluateLog = &evaluateConditionalDensityLog; }
  else { d->evaluateLog = 0; }
  d->destroy = &destroyConditionalDensity;
  d->setTransitionPoint = 0;
  d->rnd = 0;
  
  return d;
}

ConditionalDensityParameters* initializeConditionalDensityParameters(ConditionalDensityParameters* cdp, Density* density, double* fixed, bool* positions, bool *parameters, bool normalize){
  if(cdp==0) cdp = malloc(sizeof(ConditionalDensityParameters));

  cdp->normalize = normalize;
  cdp->density = density;
  cdp->fixed = fixed;
  cdp->positions = positions;
  cdp->parameters = parameters;
  
  cdp->destroyArrays = true;
  cdp->destroyDensity = true;

  return cdp;
}

double evaluateConditionalDensity(Density* this, const double* start, const double* end){
  ConditionalDensityParameters* cdp = (ConditionalDensityParameters*)this->parameters;
  
  for(int i=0; i<cdp->density->dim; i++){
    if(cdp->positions[i] == false){
      cdp->fixed[i] = *start;
      start++;
    }
  }
  
  return cdp->density->evaluate(cdp->density,cdp->fixed,cdp->fixed+cdp->density->dim);
}

double evaluateConditionalDensityLog(Density* this, const double* start, const double* end){
  ConditionalDensityParameters* cdp = (ConditionalDensityParameters*)this->parameters;
    
  for(int i=0; i<cdp->density->dim; i++){
    if(cdp->positions[i] == false){
      cdp->fixed[i] = *start;
      start++;
    }
  }
  
  return cdp->density->evaluateLog(cdp->density,cdp->fixed,cdp->fixed+cdp->density->dim);
}

void destroyConditionalDensity(Density* this){
  ConditionalDensityParameters* cdp = (ConditionalDensityParameters*)this->parameters;
  if(cdp->destroyArrays){
    free(cdp->fixed);
    free(cdp->positions);
    free(cdp->parameters);
  }
  if(cdp->destroyDensity){
    cdp->density->destroy(cdp->density);
    free(cdp->density);
  }
  free(cdp);
}
