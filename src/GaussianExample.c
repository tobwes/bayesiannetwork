#include"GaussianExample.h"

void GaussianExample(){
  int dim = 3;
  double* mu = malloc(sizeof(double)*3);
  mu[0]=1;
  mu[1]=0;
  mu[2]=3;
  double* sigma = malloc(sizeof(double)*9);
  sigma[0]=6;
  sigma[1]=5;
  sigma[2]=5;
  sigma[3]=5;
  sigma[4]=5;
  sigma[5]=4;
  sigma[6]=5;
  sigma[7]=4;
  sigma[8]=5;

  Density* gaussian = initializeGaussian(0,dim,mu,sigma);

  ProbabilitySpace* ps = initializeProbabilitySpace(0,3,gaussian);

  int n=10000;
  const double* rv[n];
  for(int i=0; i<n; i++){
    rv[i] = malloc(sizeof(double)*3);
    rv[i] = getRand(ps,true);
    printf("RV[%i]: (%f,%f,%f)\n",i,rv[i][0], rv[i][1], rv[i][2]);
  }

  double mean[dim];
  mean[0] = mean[1] = mean[2] = 0;
  for(int i=0; i<n; i++){
    mean[0] += rv[i][0];
    mean[1] += rv[i][1];
    mean[2] += rv[i][2];
  }
  mean[0] /= n;
  mean[1] /= n;
  mean[2] /= n;

  double var[dim*dim], tmp[dim];
  for(int i=0; i<dim*dim; i++){
    if(i<dim)
      tmp[i] = 0;
    var[i] = 0;
  }
  for(int i=0; i<n; i++){
    for(int j=0; j<dim; j++){
      tmp[j] = (rv[i][j] - mean[j]);
    }
    for(int j=0; j<dim*dim; j++){
      var[j] += tmp[j%3]*tmp[j/3];
    }
  }
  for(int j=0;j<dim*dim;j++){
    var[j] = var[j]/(n-1);
  }
  
  printf("mean: (%f,%f,%f)\n",mean[0],mean[1],mean[2]);
  printf("variance: (%f,%f,%f,%f,%f,%f,%f,%f,%f)\n",var[0],var[1],var[2],var[3],var[4],var[5],var[6],var[7],var[8]);  
}
