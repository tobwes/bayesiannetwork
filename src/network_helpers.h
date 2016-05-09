#ifndef NETWORK_HELPERS_H
#define NETWORK_HELPERS_H

#include"defines.h"
#include"constants.h"
#include"structures.h"
#include"includes.h"
#include"fact.h"
#include"binom.h"
#include"activation.h"
#include"helpers.h"
#include"ProbabilitySpace.h"

extern int NNUM;
extern network NETWORK;
extern neuron NEURON[];
extern double X[MAX_TRAINING_POINTS][MAX_NUM_NEURONS];
extern double Y[MAX_TRAINING_POINTS][MAX_NUM_NEURONS];
extern double OUTPUT_X[MAX_NUM_POINTS][MAX_NUM_NEURONS];
extern char OUTPUT_FILENAME[];
extern int NUMBER_OF_POINTS;
extern int OPTIMIZATION_METHOD;
extern ProbabilitySpace PROBABILITYSPACE;

void getInputVectorAsString(network *NETWORK, neuron *NEURONS, char* res);
void getWeightVectorAsString(network *NETWORK, neuron *NEURONS, char* res);
void getOutputVectorAsString(network *NETWORK, neuron *NEURONS, char* res);

void getInputVector(network *NETWORK, neuron *NEURONS, double* res);
//void getOutputVector(network *NETWORK, neuron *NEURONS, double* res);
//void getWeightsVector(network *NETWORK, neuron *NEURONS, double* res);

void saveInput(double* start); // save input state
void saveWeights(double* start);

void assignInput(const double* start); // assign input values
void assignTrainingInput(int n); // assign input values X[n]
void assignWeights(const double* start);

void feedforwardSpecial(const double *start, const double *weights, double *res);
void feedConnections(neuron *n);

double discriminant(neuron *n);
double linearDiscriminant(neuron *n);
double legendreDiscriminant(neuron *n);
double laguerreDiscriminant(neuron *n);
double fourierDiscriminant(neuron *n);

int getNetworkOutputDimensions(network *n);
int getNetworkInputDimensions(network *n);
int getNetworkWeightDimensions(network *n);

void getNetworkInput(network *net, int n, double *res);
void getTrainingInput(network *net, int n, double *res);
void getTrainingOutput(network *net, int n, double *res);

void evaluateNetworkInput(network *net);
void evaluateNetworkBayesian(network *net, double *input, NetworkStatistics *result);
void outputRV(const double* start, const double* end, double* res);

NetworkStatistics* setupNetworkStatistics(NetworkStatistics* ns, int dim);
void destroyNetworkStatistics(NetworkStatistics* ns);
#endif
