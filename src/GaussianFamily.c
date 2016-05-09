#include"GaussianFamily.h"

GaussianVarianceMatrix* initializeGaussianVarianceMatrix(GaussianVarianceMatrix* gvm, int dim, double *sigma, int factorized){

  if(gvm == 0) gvm = malloc(sizeof(GaussianVarianceMatrix));
  
  gvm->dim = dim;
  gvm->det = 1;
  
  if(!factorized){
    char equed; // equalizer info
    double rcond; // condition of matrix A

    gvm->sigma = sigma;
    gvm->S = malloc(sizeof(double)*dim);
    gvm->AF = malloc(sizeof(double)*dim*dim);
    
    LAPACKE_dposvx(LAPACK_COL_MAJOR,'E','U',gvm->dim,0,gvm->sigma,gvm->dim,gvm->AF,gvm->dim,&equed,gvm->S,0,dim,0,dim,&rcond,0,0);

    if(equed == 'N'){
      //printf("No equilibration used.\n");
      
      free(gvm->S);
      gvm->S = 0;

      for(int i = 0; i<gvm->dim; ++i){
	gvm->det*=gvm->AF[i*(dim+1)]*gvm->AF[i*(dim+1)];
      }
    }
    else {// if(equed=='Y')
      //printf("Equilibration used.\n");
      for(int i = 0; i<gvm->dim; ++i)
	gvm->det*=gvm->AF[i*(dim+1)]*gvm->AF[i*(dim+1)]/gvm->S[i]/gvm->S[i];
    }
  }else{ // if(factorized)
    gvm->sigma = 0;
    gvm->S = 0;
    gvm->AF = sigma;

    for(int i = 0; i<gvm->dim; ++i){
      gvm->det*=gvm->AF[i*(dim+1)]*gvm->AF[i*(dim+1)];
    }
  }

  return gvm;
}

GaussianParameters* initializeGaussianParameters(GaussianParameters* gp, int dim, const double* mu, double* sigma){
  if(gp==0) gp=malloc(sizeof(GaussianParameters));

  gp->dim = dim;
  
  gp->mu = mu;
  gp->free_mu = true;
  
  gp->sigma = initializeGaussianVarianceMatrix(0,dim,sigma,0);
  gp->destroy_sigma = true;

  return gp;
}

Density* initializeGaussian(Density* d, int dim, const double* mu, double* sigma){

  if(d==0) d=malloc(sizeof(Density));

  GaussianParameters* gp = initializeGaussianParameters(0,dim,mu,sigma);
  d->parameters = gp;
  
  d->dim = dim;
  d->normalizing_factor = 1; // not used
  d->log_normalizing_factor = -1./2.*log(pow(2*Pi,dim)*gp->sigma->det); // used in evaluateGaussianLog
  
  d->evaluate = &evaluateGaussian;
  d->evaluateLog = &evaluateGaussianLog;
  d->destroy = &destroyGaussian;
  d->rnd = &gaussianRnd;
  d->setTransitionPoint = &gaussianTransition;
  
  return d;
}

void solveGaussianLinearEquation(GaussianVarianceMatrix* gvm, double* rhs, int NRHS, double* x){
  char equi;
  if(gvm->S==0)
    equi = 'N';
  else
    equi = 'Y';
  double rcond;
  double ferr=1, berr=1;
  LAPACKE_dposvx(LAPACK_COL_MAJOR,'F','U',gvm->dim,NRHS,gvm->sigma,gvm->dim,gvm->AF,gvm->dim,&equi,gvm->S,rhs,gvm->dim,x,gvm->dim,&rcond,&ferr,&berr);
  //printf("ferr: %f, berr: %f\n",ferr,berr);
}


double evaluateGaussian(Density* this, const double* first, const double* last){
  Density* density = (Density*)this;
  GaussianParameters* par = (GaussianParameters*) density->parameters;
  
  double c = 1/sqrt(pow(2*Pi,par->dim)*par->sigma->det);
  
  double e = 0;
  double *x = malloc(sizeof(double)*par->dim);
  double *b = malloc(sizeof(double)*par->dim);
  
  for(int i=0;i<par->dim;i++)
    b[i] = first[i]-par->mu[i];
      
  solveGaussianLinearEquation(par->sigma,b,1,x);
  
  for(int i=0; i<par->dim; i++){
    if(par->sigma->S != 0){
      e += x[i]*b[i]/par->sigma->S[i];
    }else{
      e += x[i]*b[i];
    }
  }
  e *= -1./2.;

  return c*exp(e);
}

double evaluateGaussianLog(Density* this, const double* first, const double* last){
  GaussianParameters* par = (GaussianParameters*) this->parameters;
  
  double e = 0;
  double *x = malloc(sizeof(double)*par->dim);
  double *b = malloc(sizeof(double)*par->dim);
  
  for(int i=0;i<par->dim;i++)
    b[i] = first[i]-par->mu[i];

  solveGaussianLinearEquation(par->sigma,b,1,x); // ?!? what happens here ?!?
  
  for(int i=0; i<par->dim; i++){
    if(par->sigma->S != 0){
      e += x[i]*b[i]/par->sigma->S[i];
    }else{
      e += x[i]*b[i];
    }
  }
  e *= -1./2.;
  //double p = e+this->log_normalizing_factor;
  return e+this->log_normalizing_factor;
}

void destroyGaussianVarianceMatrix(GaussianVarianceMatrix* gvm){
  if(gvm->sigma){
    free(gvm->sigma);
    gvm->sigma = 0;
  }
  if(gvm->S){
    free(gvm->S);
    gvm->S = 0;
  }
  if(gvm->AF){
    free(gvm->AF);
    gvm->AF = 0;
  }
}

void destroyGaussian(Density* d){
  if(d->parameters){
    GaussianParameters* gp = (GaussianParameters*)d->parameters;
    destroyGaussianParameters(gp);
    free(gp);
    d->parameters = 0;
  }
}

void destroyGaussianParameters(GaussianParameters* gp){
  if(gp->destroy_sigma && gp->sigma){
    destroyGaussianVarianceMatrix(gp->sigma);
    free(gp->sigma);
    gp->sigma = 0;
  }
  if(gp->free_mu && gp->mu){
    free((double*)gp->mu);
    gp->mu = 0;
  }
}

/** maybe simplify API: remove argument const double* last **/
void gaussianTransition(Density* d, const double* first, const double* last){
  if( last - first + 1 != d->dim ){
    printf("*** Critical error: The transition point does not have the same dimension as this gaussian density.");
    exit(0);
  }

  GaussianParameters* gp = (GaussianParameters*)d->parameters;
  if(gp->free_mu && gp->mu){
    free((double*)gp->mu);
  }
  gp->mu = first;
}

/** Use Cholesky decomposition:
    sigma = t(US^{-1})US^{-1}
    -> X = mu + t(U)*Y
    -> EX = mu, VX = V(t(U)Y) = t(U)U **/
double* gaussianRnd(Density* gf, double* rv){
  GaussianParameters* gp = (GaussianParameters*)gf->parameters;
  double* y = malloc(sizeof(double)*gp->dim); // stores gp->dim independent gaussian RVs

  int r=0;
  if(gp->dim%2){ r++; }
  r+=gp->dim/2;
  
  
  double u1, u2;
  for(int i=0; i<r; i++){
    u1=rand()/(double)RAND_MAX;
    u2=rand()/(double)RAND_MAX;
    y[2*i]=sqrt(-2.0*log(u1))*cos(2.0*Pi*u2);
    if(i!=r-1 || (gp->dim%2==0))
      y[2*i+1]=sqrt(-2.0*log(u1))*sin(2.0*Pi*u2);
  }

  for(int i=0; i<gp->dim;i++){
    rv[i] = gp->mu[i];
    for(int j=0;j<i+1;j++){
      if(gp->sigma->S == 0)
	rv[i] += gp->sigma->AF[gp->dim*i+j]*y[j];
      else
	rv[i] += gp->sigma->AF[gp->dim*i+j]/gp->sigma->S[i]*y[j];
    }
  }
  
  free(y);
  return rv;
}
