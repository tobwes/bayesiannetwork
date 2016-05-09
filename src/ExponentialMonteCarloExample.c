#include"ExponentialMonteCarloExample.h"

void ExponentialMonteCarloExample(){
  double lambda = 2.3;
  Density* exponential = initializeExponential(0,lambda);

  double mu=0, sigma=2; //first try parameters
  Density* gaussian_prototype = initializeGaussian(0,1,&mu, &sigma);

  ProbabilitySpace* ps = initializeProbabilitySpace(0,1,exponential);
  setupMonteCarloEngine(ps,gaussian_prototype,METROPOLIS_SAMPLING,1000000,20);

  int n=100000;
  printf("Generating %i random variables and calculating the mean (using expectedValueFunctional)...\n",n);
  double ev_mean = expectedValueFunctional(ps, &id,n,true);
  printf("expectedValueFunctional mean:%f\n",ev_mean);

  /** Getting a pointer to the previously generated rv's **/
  const double *rv = malloc(sizeof(double)*n);
  rv = getIID(ps,n,false);

  /** Calculate the mean manually **/
  double mean = 0;
  for(int i=0; i<n; i++){
    mean+=rv[i];
  }
  mean /= n;

  /** Calculate the variance manually **/
  double var = 0;
  for(int i=0;i<n;i++){
    var += (rv[i]-mean)*(rv[i]-mean);
  }
  var /= (n-1);

  printf("Manually calculated mean and variance:\n");
  printf("mean: %f, variance: %f\n",mean,var);
}
