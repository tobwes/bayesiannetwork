/* error.c -- This belongs to gneural_network

   gneural_network is the GNU package which implements a programmable neural network.

   Copyright (C) 2016 Jean Michel Sellier
   <jeanmichel.sellier@gmail.com>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

// compute the specified error of a (training) neural network

#include "error.h"

double error(int type){
  register int n;
  double err = 0;
  
  switch(type){
    // [Linf norm] <- This is L1, isn't it?
  case LINF: // <--- please adjust
    for(n=0;n<NDATA;n++){
      // assign training input
      assignTrainingInput(n);
      
      feedforward();
      
      // compute the [L_inf] error comparing with training output <- L1?
      err+=computeLINFError(n);
    }
    return(err);
    break;
  case L2SQ: // L2 norm squared
    for(n=0;n<NDATA;n++){
      assignTrainingInput(n);
      
      feedforward();
      
      // compute the L2 error comparing with the training output for this input
      err += computeL2SQError(n);
    }
    return(err);
    break;
    
  case L2: //L2 norm
    return sqrt(error(L2SQ));
  default:
    return(0.);
    break;
  }
}

// error for one specific input
double specificError(int type, const double* start, const double* target, const double* weights, bool savestate){
  double *savedinput;
  double *savedweights;
  if(savestate){
    savedinput = malloc(MAX_IN*MAX_NUM_NEURONS*sizeof(double));
    savedweights = malloc(sizeof(double)*MAX_NUM_NEURONS*MAX_IN);
    saveInput(savedinput);
    saveWeights(savedweights);
  }

  double err=computeSpecificError(type,start,target,weights);

  if(savestate){
    assignInput(savedinput);
    assignWeights(savedweights);
    free(savedinput);
    free(savedweights);
    feedforward();
  }
  return err;
}

double computeSpecificError(int type, const double* start, const double* target, const double* weights){
  
  double output[NETWORK.num_of_neurons[NETWORK.num_of_layers-1]];
  feedforwardSpecial(start,weights,output);
  
  switch(type){
  case LINF:
    return computeSpecificLINFError(output, target);
  case L2SQ:
    return computeSpecificL2SQError(output, target);
  case L2:
    return sqrt(computeSpecificL2SQError(output,target));
  default:
    return 0;
  }
}


//only for training data X[n], Y[n]
double computeL2SQError(int n){
  double res = 0, y;
  int dimOut = getNetworkOutputDimensions(&NETWORK);
  double output[dimOut];
  getTrainingOutput(&NETWORK,n,output);
  for(int j=0; j<dimOut; ++j){
    y=NEURON[NETWORK.neuron_id[NETWORK.num_of_layers-1][j]].output;
    res+=pow(y-output[j],2);
  }
  return res;
}

//only for training data X[n], Y[n]
double computeLINFError(int n){
  double res = 0, y;
  int dimOut = getNetworkOutputDimensions(&NETWORK);
  double output[dimOut];
  getTrainingOutput(&NETWORK,n,output);
  for(int j=0; j<dimOut; ++j){
    y=NEURON[NETWORK.neuron_id[NETWORK.num_of_layers-1][j]].output;
    res+=fabs(y-output[j]);
  }
  return res;
}

// only for specific input/target combination
double computeSpecificL2SQError(double* output, const double* target){
  double res = 0;
  for(int j=0; j<NETWORK.num_of_neurons[NETWORK.num_of_layers-1]; ++j){
    res+=pow(output[j]-target[j],2);
  }
  return res;
}

// only for specific input/target combination
double computeSpecificLINFError(double* output, const double* target){
  double res = 0;
  for(int j=0;j<NETWORK.num_of_neurons[NETWORK.num_of_layers-1]; ++j){
    res+=fabs(output[j]-target[j]);
  }
  return res;
}
