#include"network_helpers.h"

void getInputVectorAsString(network *NETWORK, neuron *NEURONS, char* res){
  res[0] = '(';
  res[1] = '\0';
  int neuron_id;
  char sdouble[12];
  for(int neuron=0; neuron<NETWORK->num_of_neurons[0]; ++neuron){
     neuron_id=NETWORK->neuron_id[0][neuron];
     for(int input=0;input<NEURONS[neuron_id].nw;++input){
       sprintf(sdouble,"%.5g",NEURONS[neuron_id].x[input]);
       strcat(res,sdouble);
       if(input<NEURONS[neuron_id].nw-1){ strcat(res,", "); }
     }
     if(neuron<NETWORK->num_of_neurons[0]-1){ strcat(res,", "); }
   }
  strcat(res,")");
}

void getInputVector(network *NETWORK, neuron *NEURONS, double* res){
  int neuron_id;
  for(int neuron=0; neuron<NETWORK->num_of_neurons[0]; ++neuron){
    neuron_id = NETWORK->neuron_id[0][neuron];
    for(int input=0; input<NEURONS[neuron_id].nw; ++input){
      *res = NEURONS[neuron_id].x[input];
      ++res;
    }
  }
}

void getOutputVectorAsString(network *NETWORK, neuron *NEURONS, char* res){
  res[0] = '(';
  res[1] = '\0';
  int neuron_id;
  char sdouble[12];
  for(int neuron=0; neuron<NETWORK->num_of_neurons[NETWORK->num_of_layers-1]; ++neuron){
     neuron_id=NETWORK->neuron_id[NETWORK->num_of_layers-1][neuron];
     sprintf(sdouble,"%.5g",NEURONS[neuron_id].output);
     strcat(res,sdouble);
     if(neuron<NETWORK->num_of_neurons[NETWORK->num_of_layers-1]-1){ strcat(res,", "); }
  }
  strcat(res,")");
}

void getWeightVectorAsString(network *NETWORK, neuron *NEURONS, char* res){
  res[0] = '(';
  res[1] = '\0';
  int neuron_id;
  char sdouble[12];
  for(int layer=1; layer<NETWORK->num_of_layers; ++layer){
    for(int neuron=0; neuron<NETWORK->num_of_neurons[layer]; ++neuron){
      neuron_id=NETWORK->neuron_id[layer][neuron];
      for(int input=0; input<NEURONS[neuron_id].nw; ++input){
	sprintf(sdouble,"%.5g",NEURONS[neuron_id].w[input]);
	strcat(res,sdouble);
	if(input<NEURONS[neuron_id].nw-1){ strcat(res,", "); }
      }
      if(layer<NETWORK->num_of_layers-1 || neuron<NETWORK->num_of_neurons[layer]-1){ strcat(res,", "); }
    }
  }
  strcat(res,")");
}

void saveInput(double* start){
  int neuron_id;
  for(int local_id=0; local_id<NETWORK.num_of_neurons[0]; ++local_id){
    neuron_id=NETWORK.neuron_id[0][local_id];
    *start = NEURON[neuron_id].x[0];
    ++start;
  }
}

void saveWeights(double* start){
  int neuron_id;
  for(int layer=1; layer<NETWORK.num_of_layers; ++layer){
    for(int local_id=0; local_id<NETWORK.num_of_neurons[layer]; ++local_id){
      neuron_id=NETWORK.neuron_id[layer][local_id];
      for(int j=0; j<NEURON[neuron_id].nw; ++j){
	*start = NEURON[neuron_id].w[j];
	++start;
      }
    }
  }
}

void assignInput(const double* start){
  int neuron_id;
  for(int local_id=0; local_id<NETWORK.num_of_neurons[0]; ++local_id){
    neuron_id=NETWORK.neuron_id[0][local_id];
    NEURON[neuron_id].x[0]=NEURON[neuron_id].output=*start;
    ++start;
  }
}

void assignTrainingInput(int n){
  int neuron_id;
  for(int local_id=0; local_id<NETWORK.num_of_neurons[0]; ++local_id){
    neuron_id=NETWORK.neuron_id[0][local_id];
    NEURON[neuron_id].x[0]=NEURON[neuron_id].output=X[n][neuron_id];
  }
}

void assignWeights(const double* start){
  for(int layer=1; layer<NETWORK.num_of_layers; ++layer){
    for(int local_id=0; local_id<NETWORK.num_of_neurons[layer]; ++local_id){
      int neuron_id = NETWORK.neuron_id[layer][local_id];
      for(int connection=0; connection<NEURON[neuron_id].nw; ++connection){
	NEURON[neuron_id].w[connection]=*start;
	++start;
      }
    }
  }
}

void feedforwardSpecial(const double *input, const double *weights, double *res){

  if(weights!=0){ assignWeights(weights); }
  if(input!=0){ assignInput(input); }
      
  for(int layer=0; layer<NETWORK.num_of_layers; layer++){
    for(int local_id=0; local_id<NETWORK.num_of_neurons[layer]; ++local_id){
      int neuron_id = NETWORK.neuron_id[layer][local_id];
      
      if(layer>0){
	feedConnections(&(NEURON[neuron_id]));
      }
      
      double discr = discriminant(&(NEURON[neuron_id]));
      
      NEURON[neuron_id].output = activation(NEURON[neuron_id].activation, discr);
      
      if(res && layer == NETWORK.num_of_layers-1){
	res[local_id]=NEURON[neuron_id].output;
      }
    }  
  }
}

double discriminant(neuron *n){
      switch(n->discriminant){
      case LINEAR:
	return linearDiscriminant(n);
	break;
      case LEGENDRE:
	return legendreDiscriminant(n);
	break;
      case LAGUERRE:
	return laguerreDiscriminant(n);
	break;
      case FOURIER:
	return fourierDiscriminant(n);
	break;
      default:
	return 0;
	break;
      }
}

double linearDiscriminant(neuron *n){
  double x=0;
  for(int i=0; i<n->nw; ++i){
    x+=n->x[i]*n->w[i]; // linear product between w[] and x[]
  }
  return x;
}

double legendreDiscriminant(neuron *n){
  double disc=0;
  double a,tmp;
  for(int i=0; i<n->nw; ++i){
    a=pow(2,i);
    tmp=0.;
    for(int j=0; j<=i; ++j){
      tmp+=pow(n->x[i],j)*binom(i,j)*binom((i+j-1)/2,j);
    }
    tmp*=a*n->w[i];
    disc+=tmp;
  }
  return disc;
}

double laguerreDiscriminant(neuron *n){
  double disc=0;
  double tmp;
  for(int i=0; i<n->nw; i++){
      tmp=0.;
      for(int j=0;j<=i;j++){
	tmp+=binom(i,j)*pow(n->x[i],j)*pow(-1,j)/fact(j);
      }
      tmp*=n->w[i];
      disc+=tmp;
  }
  return disc;
}

double fourierDiscriminant(neuron *n){
  double disc=0;
  double tmp;
  for(int i=0; i<n->nw;i++){
    tmp=0.;
    for(int j=0;j<=i;j++){
      tmp+=sin(2.*j*Pi*n->x[i]);
    }
    tmp*=n->w[i];
    disc+=tmp;
  }
  return disc;
}

void feedConnections(neuron *n){
  for(int connection=0; connection<n->nw; ++connection){
    n->x[connection]=NEURON[n->connection[connection]].output;
  }
}

int getNetworkOutputDimensions(network *n){
  int dim = n->num_of_neurons[n->num_of_layers-1];
  return dim;
}

int getNetworkInputDimensions(network *n){
  int dim = n->num_of_neurons[0];
  return dim;
}

int getNetworkWeightDimensions(network *n){
  int dim = 0;
  for(int layer=1; layer<NETWORK.num_of_layers; ++layer){
    for(int local_id=0; local_id<NETWORK.num_of_neurons[layer]; ++local_id){
      int neuron_id = NETWORK.neuron_id[layer][local_id];
      dim += NEURON[neuron_id].nw;
    }
  }
  return dim;
}

void evaluateNetworkInput(network *net){
  // open file
  FILE *fp=fopen(OUTPUT_FILENAME,"w");

  //int dimWeights = getNetworkWeightDimensions(net);
  int dimOut = getNetworkOutputDimensions(net);
  int dimIN = getNetworkInputDimensions(net);
  //double weights[dimWeights];
  char inputstring[dimIN*2+2];
  
  for(int n=0; n<NUMBER_OF_POINTS; ++n){
    NetworkStatistics *ns = setupNetworkStatistics(0,dimOut);
    ns->copy_rvs = true;
    ns->output_rvs = malloc(sizeof(double)*dimOut*10000);
    
    // assign input
    getNetworkInput(net,n,net->input);

    // print input to terminal
    vectorToString(net->input,dimIN,inputstring);
    printf("Evaluating input x=%s\n",inputstring);

    // print input to file
    for(int i=0; i<dimIN; ++i){
      fprintf(fp,"%g,",net->input[i]);
    }
    
    // feedforward
    if(OPTIMIZATION_METHOD == BAYESIAN){
      evaluateNetworkBayesian(net,net->input,ns);
    }else{
      feedforwardSpecial(net->input,0,ns->mean);
    }
    
    // print mean output
    for(int i=0; i<dimOut; ++i){
      if(i<dimOut-1){
	fprintf(fp,"%g,",ns->mean[i]);
      }else{
	fprintf(fp,"%g",ns->mean[i]);
      }
    }

    if(OPTIMIZATION_METHOD == BAYESIAN){
      fprintf(fp,",");
      //print variance output
      for(int i=0; i<dimOut*dimOut; ++i){
	if(i<dimOut-1 || dimOut ==1){
	  fprintf(fp,"%g,",ns->variance[i]);
	}else{
	  fprintf(fp,"%g\n",ns->variance[i]);
	}
      }
      if(dimOut == 1){
	fprintf(fp,"%g,%g,%g,%g,%g,%g\n",ns->sd,ns->lower_quantile,ns->lower_quartile,ns->median,ns->upper_quartile,ns->upper_quantile);
      }
    }else{
      fprintf(fp,"\n");
    }
    destroyNetworkStatistics(ns);
    free(ns);
  }
}

NetworkStatistics* setupNetworkStatistics(NetworkStatistics* ns, int dim){
  if(ns == 0) { ns = malloc(sizeof(NetworkStatistics)); }

  ns->mean = malloc(sizeof(double)*dim);
  ns->variance = malloc(sizeof(double)*dim*dim);

  ns->output_rvs = 0;
  ns->copy_rvs = false;
  return ns;
}

void destroyNetworkStatistics(NetworkStatistics* ns){
  if(ns->mean != 0) { free(ns->mean); }
  if(ns->variance != 0) { free(ns->variance); }
  if(ns->copy_rvs && ns->output_rvs != 0) { free((double*)ns->output_rvs); }
}

void evaluateNetworkBayesian(network *net, double *input, NetworkStatistics *result){
  int dimOut = getNetworkOutputDimensions(net);
  int dimWeights = getNetworkWeightDimensions(net);
  int nrvs = 10000;
  const double* weight_rvs;
  weight_rvs = getIID(&PROBABILITYSPACE, nrvs, false, 0);
  for(int i=0; i<nrvs; ++i){
    feedforwardSpecial(net->input,weight_rvs+i*dimWeights,(double*)(result->output_rvs+i*dimOut));
  }

  mean(result->output_rvs,dimOut,nrvs,result->mean);
  variance(result->output_rvs,dimOut,nrvs,result->mean,result->variance);
  if(dimOut == 1){
    double q[5];
    q[0] = 0.2;
    q[1] = 0.25;
    q[2] = 0.5;
    q[3] = 0.75;
    q[4] = 0.98;
    double q_res[5];
    quantile(result->output_rvs,nrvs,q,5,q_res);
    result->lower_quantile = q_res[0];
    result->lower_quartile = q_res[1];
    result->median = q_res[2];
    result->upper_quartile = q_res[3];
    result->upper_quantile = q_res[4];

    result->sd = sqrt(*(result->variance));
  }
  
  for(int j=0; j<dimOut; ++j){
    NEURON[net->neuron_id[net->num_of_layers-1][j]].output = result->mean[j];
  }
}

void outputRV(const double* start, const double* end, double* res){
  feedforwardSpecial(NETWORK.input, start, res);
}

void getNetworkInput(network *net, int n, double *res){
  for(int i=0; i<getNetworkInputDimensions(net); ++i){
    res[i] = OUTPUT_X[n][i];
  }
}

void getTrainingInput(network *net, int n, double *res){
  for(int i=0; i<getNetworkInputDimensions(net); ++i){
    res[i] = X[n][i];
  }
}

void getTrainingOutput(network *net, int n, double *res){
  for(int i=0; i<getNetworkOutputDimensions(net); ++i){
    res[i] = Y[n][i];
  }
}
