#include "bayesian.h"

void setup_probabilityspace(int verbosity){

  /** Determining dimensions **/
  int dimWeights = getNetworkWeightDimensions(&NETWORK);  // weight vector dimension; read dim from network
  int dimX = getNetworkInputDimensions(&NETWORK);  // input dimension
  int dimY = getNetworkOutputDimensions(&NETWORK); // output dimension == number of neurons in the last layer
  int dimXYW = dimX + dimY + dimWeights; // total dimension = input dim + output dim + weights dim

  /** Creating prior density **/
  double *mu_prior = rep(0,dimWeights,0), *sigma_prior=diag(100,dimWeights,dimWeights,0); //prior estimates
  Density* weight_prior = initializeGaussian(0,dimWeights,mu_prior,sigma_prior); // need not be a gaussian

  /** Error density is the (common) density of (x,y,theta) **/
  double beta = 1;
  Density* error_density = initializeErrorDensity(0,ERROR_TYPE,dimX,dimY,dimWeights,beta);
  ErrorDensityParameters *edp = (ErrorDensityParameters*) error_density->parameters;

  
  Density* data_conditional = 0; // pointer to density error_density, with fixed parameters (x,y)
  Density* data_conditionals_product = 0; // pointer to the complete product density

  bool positions[dimXYW]; // an array with the fixed positions (will be copied for each data_conditional density)
  for(int i=0; i<dimX+dimY; ++i){ positions[i] = true; } // 
  for(int i=0; i<dimWeights; ++i){ positions[dimX+dimY+i] = false; }

  for(int i=0; i<NDATA; ++i){
    double *fixed;
    bool* positions_cpy;
    
    fixed = malloc(sizeof(double)*(dimX+dimY+dimWeights));

    getTrainingInput(&NETWORK,i,fixed);
    getTrainingOutput(&NETWORK,i,fixed+dimX);
    
    positions_cpy = malloc(sizeof(bool)*dimXYW);
    memcpy(positions_cpy,positions,sizeof(bool)*dimXYW);

    if(i==0){
      data_conditional = data_conditionals_product = initializeConditionalDensity(0, error_density, fixed, positions_cpy, 0, false);
    }else{
      data_conditional = initializeConditionalDensity(0, error_density, fixed, positions_cpy, 0, false);
      ConditionalDensityParameters* cdp = (ConditionalDensityParameters*)data_conditional->parameters;
      cdp->destroyDensity = false;

      data_conditionals_product = initializeProductDensity(0,data_conditionals_product,data_conditional);
    }
  }
  
  /** Creating posterior density: product of data_conditionals and prior **/
  Density* weight_posterior = initializeProductDensity(0,weight_prior,data_conditionals_product);

  /** setting up probability space **/
  initializeProbabilitySpace(&PROBABILITYSPACE,dimWeights,weight_posterior);
  for(int layer=1; layer<NETWORK.num_of_layers; ++layer){
    for(int local_id=0; local_id<NETWORK.num_of_neurons[layer]; ++local_id){
      int neuron_id = NETWORK.neuron_id[layer][local_id];
      registerRVs(&PROBABILITYSPACE,NEURON[neuron_id].w,NEURON[neuron_id].nw);
    }
  }

  /** Setting up Monte Carlo Enginge **/
  double *mu_proto = rep(0,dimWeights,0), // needs to be free'd after burnin
    *sigma_proto=diag(2,dimWeights,dimWeights,0); // will be free'd by gaussian_prototype
  Density* gaussian_prototype = initializeGaussian(0,dimWeights,mu_proto,sigma_proto);
  GaussianParameters* gp = (GaussianParameters*) gaussian_prototype->parameters;
  gp->free_mu = false; // mu cannot be free'd, because it will be set to the random transition points, which are cached in the Monte Carlo engine
  
  setupMonteCarloEngine(&PROBABILITYSPACE,gaussian_prototype,METROPOLIS_SAMPLING,10000,1);
  new_burnin(&PROBABILITYSPACE,true);
  free(mu_proto);
  mu_proto = 0;

  double err = 1.e8; // just a huge error
  int runs = MAX_OPTIMIZATION_RUNS;
  printf("Calculating network error...\n");
  while( (err = calculate_training_error()) > ERROR_BOUND && runs > 0 ){
    double new_beta = (1+err)*edp->beta;
    printf("The network error is %g. This is too much.\n\
Adjusting beta. New beta: %f\n\
Calculating network error...\n",err,new_beta);
    set_beta(error_density,new_beta);
  }

  if(err > ERROR_BOUND){
    printf("*** Critical error: could not optimize the network weights.");
    exit(0);
  }

  printf("The network error is %g. This is good enough. The network is now set up.\n",err);
}

double calculate_training_error(){  
  double err;

  new_burnin(&PROBABILITYSPACE,false);
  expectedValue(&PROBABILITYSPACE,&errorFunction,1,10000,&err);
  return err;
}

void errorFunction(double* res){
  *res = error(ERROR_TYPE);
}
