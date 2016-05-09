#include"ProductDensity.h"

Density* initializeProductDensity(Density* d, Density* d1, Density* d2){
  if(d==0){ d=malloc(sizeof(Density)); }

  d->parameters = initializeProductDensityParameters(0,d1,d2);

  if( d1->dim != d2->dim ){
    printf("ProductDensity.c: The dimensions of these densities are not the same!");
    exit(0);
  }
  
  d->dim = d1->dim;
  d->normalizing_factor = 1;
  d->log_normalizing_factor = 0;
  d->evaluate = &evaluateProductDensity;
  if(d1->evaluateLog && d2->evaluateLog) { d->evaluateLog = &evaluateProductDensityLog; }
  else { d->evaluateLog = 0; }
  d->destroy = &destroyProductDensity;
  d->setTransitionPoint = 0;
  d->rnd = 0;
  
  return d;
}

ProductDensityParameters* initializeProductDensityParameters(ProductDensityParameters* pdp, Density* d1, Density* d2){
  if(pdp == 0) pdp = malloc(sizeof(ProductDensityParameters));

  pdp->d1 = d1;
  pdp->d2 = d2;
  pdp->destroyD1 = true;
  pdp->destroyD2 = true;
  
  return pdp;
}

double evaluateProductDensity(Density* this, const double* start, const double* end){
  ProductDensityParameters* pdp = (ProductDensityParameters*)this->parameters;

  if(this->evaluateLog){
    double logp =this->evaluateLog(this,start,end);
    double p=exp(logp);
    return p;
  }else{
    double p1 = pdp->d1->evaluate(pdp->d1,start,end), p2 = pdp->d2->evaluate(pdp->d2,start,end), p=p1*p2;
    return p;
  }
}

double evaluateProductDensityLog(Density* this, const double* start, const double* end){
  ProductDensityParameters* pdp = (ProductDensityParameters*)this->parameters;
  
  double p1,p2;
  p1 = pdp->d1->evaluateLog(pdp->d1,start,end);
  p2 = pdp->d2->evaluateLog(pdp->d2,start,end);
  double p = p1+p2;
  return p;
}

void destroyProductDensity(Density* this){
  ProductDensityParameters* pdp = (ProductDensityParameters*)this->parameters;
  if(pdp->destroyD1 && pdp->d1){
    pdp->d1->destroy(pdp->d1);
    free(pdp->d1);
    pdp->d2 = 0;
  }
  if(pdp->destroyD2 && pdp->d2){
    pdp->d2->destroy(pdp->d2);
    free(pdp->d2);
    pdp->d2 = 0;
  }
  free(pdp);
}
