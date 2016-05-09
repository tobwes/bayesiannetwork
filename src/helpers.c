#include "helpers.h"

void vectorToString(const double* vec, int dim, char* str){
  str[0] = '(';
  str[1] = '\0';
  char tmp[20] = "";
  for(int i=0; i<dim; ++i){
    if(i<dim-1){
      sprintf(tmp,"%g,",vec[i]);
    }
    else if(i==dim-1){
      sprintf(tmp,"%g)\n",vec[i]);
    }
    strcat(str,tmp);
  }
}

void mean(const double* rvs, int dim, int nrvs, double* m){
  for(int i=0; i<dim; ++i){
    m[i] = 0;
  }
  const double* rvs_end = rvs + dim*nrvs;
  while(rvs != rvs_end){
    for(int i = 0; i<dim; ++i){
      m[i] += rvs[i]/nrvs;
    }
    rvs += dim;
  }
}

void variance(const double* rvs, int dim, int nrvs, const double *m, double *var){
  const double *mm = m; // working copy of the mean

  // if no mean provided, calculate mean
  if(mm == 0){
    mm = malloc(sizeof(double)*dim);
    mean(rvs,dim,nrvs,(double*)mm);
  }

  // set var to 0;
  for(int i=0; i<dim; ++i){
    for(int j=0; j<dim; ++j){
      var[j+dim*i] = 0;
    }
  }

  // calculate variance
  for(const double* rvs_end = rvs+dim*nrvs; rvs<rvs_end; rvs+=dim){
    for(int i=0; i<dim; ++i){
      for(int j=0; j<dim; ++j){
	var[j+dim*i] += (rvs[j]*rvs[i]- mm[j]*mm[i])/nrvs;
      }
    }
  }
}

int double_leq(const void* p1, const void* p2){
  double d = *(double*)p1 - *(double*)p2;

  if(d>0){ return 1; }
  else if(d<0){ return -1; }
  else return 0;
}

void quantile(const double* rvs, int nrvs, double *q, int n_quantiles, double* res){
  double* rvs_copy = malloc(sizeof(double)*nrvs);
  memcpy(rvs_copy,rvs,sizeof(double)*nrvs);

  qsort(rvs_copy,nrvs,sizeof(double),&double_leq);

  for(int i=0; i<n_quantiles; ++i){
    res[i] = rvs_copy[(int)floor(q[i]*(double)nrvs)];
  }
  free(rvs_copy);
}
