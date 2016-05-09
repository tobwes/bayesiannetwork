#include"lapack_helpers.h"

double* rep(double d, int n, double* res){
  if(res == 0){ res = malloc(sizeof(double)*n); }

  for(int i=0; i<n; ++i){
    res[i] = d;
  }

  return res;
}

double* diag(double d, int m, int n, double* res){
  if(res == 0){ res = malloc(sizeof(double)*m*n); }

  for(int i=0; i<m; ++i){
    for(int j=0; j<n; ++j){
      if(i==j){
	*(res+n*i+j) = d;
      }else{
	*(res+n*i+j) = 0;
      }
    }
  }

  return res;
}

void printArray(double* array, int m, int n){
  printf("Array output:\n(");
  for(int j=0; j<n; ++j){
    for(int i=0; i<m ; ++i){
      printf("%g",*(array+n*i+j));
      if(i<m-1)
	printf(",");
    }
    if(j<n-1)
      printf(";\n");
  }
  printf(")\n\n");
}
