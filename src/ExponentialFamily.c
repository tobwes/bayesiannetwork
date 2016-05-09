#include"ExponentialFamily.h"

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
						       void (*T)(const double* first, const double* last, double* res)){
  if(ep==0) ep=malloc(sizeof(ExponentialParameters));

  ep->dim = dim;
  ep->dimTheta = dimTheta;
  ep->dimScalarProduct = dimScalarProduct;
  ep->theta = theta;
  ep->free_theta = true;
  ep->_g = _g;
  ep->_A = _A;
  ep->_eta = malloc(sizeof(double)*dim);
  for(int i=0; i<dimScalarProduct; ++i){
    ep->_eta[i] = _eta ? _eta[i] : nan("");
  }
  ep->g = g;
  ep->A = A;
  ep->eta = eta;
  ep->h = h;
  ep->T = T;

  if(g && isnan(_g) && theta){
    g(theta,theta+dimTheta, &(ep->_g));
  }
  if(A && isnan(_A) && theta){
    A(theta,theta+dimTheta, &(ep->_A));
  }
  if(eta && isnan(_eta[0]) && theta){
    eta(theta,theta+dimTheta,ep->_eta);
  }

  if( (isnan(_g) && isnan(_A)) || !_eta ){
    printf("*** Critical error: too few parameters in initializing exponential family\n");
    exit(0);
  }
  
  return ep;
}

Density* initializeExponentialFamily(Density* ef,
				     int dim,
				     int dimTheta,
				     int dimScalarProduct,
				     double *theta,
				     double _g,
				     double _A,
				     double *_eta,
				     void (*g)(const double* first, const double* last, double* res),
				     void (*A)(const double* first, const double* last, double* res),
				     void (*eta)(const double* first, const double* last, double* res),
				     void (*h)(const double* first, const double* last, double* res),
				     void (*T)(const double* first, const double* last, double* res)){
  if(ef==0) ef = malloc(sizeof(Density));
  ExponentialParameters* ep = initializeExponentialParameters(0,dim,dimTheta,dimScalarProduct,theta,_g,_A,_eta,g,A,eta,h,T);
  ef->parameters = ep;

  ef->dim = dim;
  ef->normalizing_factor = 1; // not used
  ef->log_normalizing_factor = 0;
  if(h == &one && ep->_g){ ef->log_normalizing_factor = log(ep->_g); }
  ef->evaluate = &evaluateExponentialFamily;
  if(h == &one){ ef->evaluateLog = &evaluateExponentialFamilyLog; }
  else{ ef->evaluateLog = 0; }
  ef->destroy = &destroyExponentialFamily;
  ef->setTransitionPoint = 0;

  return ef;
}
  
double evaluateExponentialFamily(Density* this, const double* first, const double* last){

  ExponentialParameters* ep = (ExponentialParameters*)(this->parameters);

  double t[ep->dim];

  ep->T(first,last,t);

  double sp = 0; // scalar product
  for(int i=0; i<ep->dimScalarProduct; ++i){
    sp += ep->_eta[i]*t[i];
  }

  double res = nan("");
  if(!isnan(ep->_A)){
    res = exp(sp+ep->_A);
  }else if(!isnan(ep->_g)){
    res = ep->_g*exp(sp);
  }
  
  double _h;
  ep->h(first,last,&_h);
  res *= _h;
  
  return res;
}

double evaluateExponentialFamilyLog(Density* this, const double* first, const double* last){
  ExponentialParameters* ep = (ExponentialParameters*)(this->parameters);

  double t[ep->dim];

  ep->T(first,last,t);

  double sp = 0; // scalar product
  for(int i=0; i<ep->dim; ++i){
    sp += ep->_eta[i]*t[i];
  }

  double res = nan("");
  if(!isnan(ep->_A)){
    res = sp+ep->_A;
  }else if(!isnan(ep->_g)){
    res = sp;
  }
    
  return res + this->log_normalizing_factor;
}

void destroyExponentialFamily(Density* ef){
  if(ef->parameters){
    destroyExponentialParameters((ExponentialParameters*)(ef->parameters));
    free(ef->parameters);
    ef->parameters = 0;
  }
}

void destroyExponentialParameters(ExponentialParameters* ep){
  if(ep->_eta){
    free(ep->_eta);
    ep->_eta = 0;
  }
  if(ep->theta && ep->free_theta){
    free(ep->theta);
    ep->theta = 0;
  }
}

void setTheta(Density* ef, double* theta){
  ExponentialParameters* ep = (ExponentialParameters*)(ef->parameters);
  ep->theta = theta;
    
  if(ep->g){
    ep->g(theta,theta+ep->dim, &(ep->_g));
    if(ep->h == &one){
      ef->log_normalizing_factor = log(ep->_g);
    }
  }
  if(ep->A){
    ep->A(theta,theta+ep->dimTheta, &(ep->_A));
  }
  if(ep->eta){
    ep->eta(theta,theta+ep->dimTheta,ep->_eta);
  }
  
}
 
Density* initializeExponential(Density* ed, double lambda){
  if(ed==0) ed = malloc(sizeof(Density));

  double _eta;
  _eta = -lambda;
  
  initializeExponentialFamily(ed,1,1,1,0,lambda,nan(""),&_eta,0,0,0,heaviside,id);

  return ed;
}

void destroyExponential(Density* e){
  destroyExponentialFamily(e);
}
