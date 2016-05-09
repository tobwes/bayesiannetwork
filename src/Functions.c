#include"Functions.h"

void id(const double* first, const double* last, double* res){
  int dim = last-first+1;
  for(int i=0; i<dim; ++i)
    res[i] = first[i];
}

void heaviside(const double* first, const double* last, double* res){
  *res = 1;
  int dim = last - first + 1;
  for(int i=0;i<dim;++i)
    *res = *res==1 && first[i]>=0 ? 1 : 0;
}

void one(const double* first, const double* last, double* res){
  *res = 1;
}
  
