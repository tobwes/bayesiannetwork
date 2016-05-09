#include"ProbabilitySpace.h"

ProbabilitySpace* initializeProbabilitySpace(ProbabilitySpace* ps, int dim, Density* d){
  if(ps==0){
    ps = malloc(sizeof(ProbabilitySpace));
  }

  ps->dim = dim;
  ps->end = ps->begin = malloc(sizeof(double*)*dim);
  ps->density = d;
  if(d!=0 && d->rnd != 0){
    ps->rv_cache_end = ps->rv_cache_start = malloc(sizeof(double)*dim*RV_CACHE_SIZE);
  }
  ps->mce = 0;
  ps->registered_rvs = 0;

  return ps;
}

void destroyProbabilitySpace(ProbabilitySpace* ps){
  for(int i=0; i<ps->dim - ps->registered_rvs; i++){
    if(ps->begin[ps->dim-1-i]) free(ps->begin[ps->dim-1-i]);
  }
  free(ps->begin);
  if(ps->density){
    ps->density->destroy(ps->density);
    free(ps->density);
  }
  if(ps->mce){
    destroyMonteCarloEngine(ps->mce);
    free(ps->mce);
  }
  if(ps->rv_cache_start){
    free(ps->rv_cache_start);
  }
}

void registerRVs(ProbabilitySpace* ps, double* begin, int n){
  // try to write at these RVs (just to make debuggin easier, as this will throw an error, when trying to register orphaned pointers)
  // the computational effort is minimal, too...
  for(int i=0; i<n; ++i){
    begin[i] = 0;
  }
  
  int cnt = n;
  while(cnt > 0 && ps->end - ps->begin < ps->dim){
    *(ps->end) = begin;
    (ps->end)++;
    begin++;
    cnt--;
  }

  if(cnt>0){
    printf("*** Critical error: trying to register too many RVs at a probability space\n");
    exit(0);
  }

  ps->registered_rvs += n;
}

void registerRV(ProbabilitySpace* ps, double* rv){
  registerRVs(ps, rv, 1);
}

void createRVs(ProbabilitySpace* ps){
  printf("Creating %i RVs...\n",(int)(ps->dim-(ps->end-ps->begin)));
  for(int i=0; i<ps->dim-(ps->end-ps->begin); i++){
    ps->begin[ps->dim-i-1] = malloc(sizeof(double));
  }
}

const double* getIID(ProbabilitySpace* ps, int n, bool generate_new, const double* res){
  if(ps->end-ps->begin < ps->dim){ createRVs(ps); }
  
  const double *rres; // internal result array
  if(ps->density->rnd == 0){
    MonteCarloEngine* mce = ps->mce;
    if(mce == 0){
      printf("*** Critical error: you need to set up a Monte Carlo engine before generating (iid) random variables\n");
      exit(0);
    }
    
    rres = getRVs(ps,n,generate_new);
  }else{
    int m=n;
    if(generate_new || (m = n-(int)(ps->rv_cache_end - ps->rv_cache_start)/ps->dim) > 0){
      if((int)(ps->rv_cache_end-ps->rv_cache_start)/ps->dim + m > RV_CACHE_SIZE){
	printf("*** Critical error: RV cache overflow\n");
	exit(0);
      }
      while(m>0){
	ps->density->rnd(ps->density, ps->rv_cache_end);
	ps->rv_cache_end += ps->dim;
	m--;
      }
    }
    rres = ps->rv_cache_end-ps->dim*n;
  }

  // saving last state
  const double* rres_end = rres + n*ps->dim;
  for(int i=0; i<ps->dim; i++){
    rres_end --;
    *(ps->begin[ps->dim-1-i]) = *rres_end;
  }

  if(res == 0){
    res = rres;
  }else{
    memcpy((double*)res,rres,sizeof(double)*ps->dim*n);
  }
  return res;
}

const double* getRand(ProbabilitySpace* ps, bool generate_new){
  return getIID(ps,1,generate_new,0);
}

void expectedValue(ProbabilitySpace* ps, void (*f)(double* res), int dimf, int nrvs, double* ev){
  int i = nrvs;
  *ev = 0;
  double *tmp = malloc(sizeof(double)*dimf);
  while( i > 0 ){
    getRand(ps,true);
    (*f)(tmp);
    for(int j=0; j<dimf; ++j){
      ev[j] += tmp[j]/(double)nrvs;
    }
    i--;
  }
  free(tmp);
}

void expectedValueFunctional(ProbabilitySpace* ps, void (*f)(const double* first, const double* last, double* res), int dimf, int nrvs, bool newRVs, double* res){
  for(int j=0; j<dimf; ++j){
    res[j] = 0;
  }
  const double *rvs = getIID(ps, nrvs, newRVs,0);
  int i = nrvs;
  double fres[dimf];
  while( i > 0 ){
    i--;
    (*f)(rvs+i*ps->dim, rvs+(i+1)*ps->dim-1, fres);
    for(int j=0; j<dimf; j++){
      res[j] += fres[j]/(double)nrvs;
    }
  }
}
void setupMonteCarloEngine(ProbabilitySpace* ps, Density* sampling_prototype, int sampling_algorithm, int burnin, int gap, ...){
  if(ps == 0){
    printf("*** Critical error: you should supply a probability space to associate the MonteCarloEngine with.\n");
    exit(0);
  }
  if(sampling_prototype->rnd == 0){
    printf("*** Critical error: cannot draw RVs from sampling density\n");
    exit(0);
  }

  if(sampling_algorithm != REJECTION_SAMPLING){
    ps->mce = initializeMonteCarloEngine(ps,ps->density,sampling_prototype,sampling_algorithm,burnin,gap,0);
  }else{ //if(sampling_algorithm == REJECTION_SAMPLING)
    va_list valist;
    va_start(valist, gap);
    va_arg(valist,ProbabilitySpace*);
    va_arg(valist,Density*);
    va_arg(valist, int);
    va_arg(valist, int);
    va_arg(valist, int);
    double rejection_sampling_factor = va_arg(valist,double);
    va_end(valist);
    initializeMonteCarloEngine(ps,ps->density,sampling_prototype,sampling_algorithm,burnin,gap,rejection_sampling_factor);
  }  
}

void setState(ProbabilitySpace* ps, double* start){
  for(int i=0; i<ps->dim; ++i){
    *(ps->begin[i]) = start[i];
  }
    
}
