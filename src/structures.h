/* rnd.h -- This belongs to gneural_network

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

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "defines.h"
#include "includes.h"

typedef struct _neuron{
 int nw;
 // connection[i] is the identification number of the neuron which
 // output is connected to the i-th input branch of the neuron
 int connection[MAX_IN];
 int activation; // type of activation function
 int discriminant; // type of discriminant function
 double x[MAX_IN]; // n inputs
 double w[MAX_IN]; // n weights
 double output;    // one output
} neuron;

typedef struct _network{
  int num_of_layers; // total number of layers
  int num_of_neurons[MAX_NUM_LAYERS]; // number of neurons per layer
  // neuron_id[i][j] is the global identification number of the j-th neuron in the i-th layer
  int neuron_id[MAX_NUM_LAYERS][MAX_NUM_NEURONS];
  double input[MAX_NUM_NEURONS*MAX_IN];
} network;

typedef struct _Density{
  /**---------------------------------------------------------------------------------------------------------------------------------/
  /---- Density specific parameters --------------------------------------------------------------------------------------------------/
  /---- They are owned (and free'd) by this density struct ---------------------------------------------------------------------------/
  /---------------------------------------------------------------------------------------------------------------------------------**/
  void* parameters;

  /** the dimension of the density **/
  int dim;

  /** A constant normalizing factor for the density and its log **/
  double normalizing_factor, log_normalizing_factor;

  /** a pointer to the specific evaluate function of the density **/
  double (*evaluate)(struct _Density* this, const double* first, const double* last);
  /** a pointer to the log of the evaluate function **/
  /** it should only be implemented if it is more efficient than log(evaluate()) **/
  double (*evaluateLog)(struct _Density* this, const double* first, const double* last);
  /** a pointer to a function drawing random variables directly from this density (NO Monte Carlo!)**/
  double* (*rnd)(struct _Density* d, double* start);
  /** a pointer to a function setting transition points, if this is a family of transition kernels **/
  void (*setTransitionPoint)(struct _Density* d, const double* start, const double* end);
  /** density-specific method to destroy this density **/
  void (*destroy)(struct _Density* d);
} Density;

typedef enum _SamplingAlgorithm{IMPORTANCE_SAMPLING,
				REJECTION_SAMPLING,
				METROPOLIS_SAMPLING,
				GIBBS_SAMPLING
} SamplingAlgorithm;
typedef struct _MonteCarloEngine{
  /** the gap between two used samples **/
  int gap;
  /** the burnin, when the engine is started **/
  int burnin;

  /** the used sampling algorithm **/
  SamplingAlgorithm sampling_algorithm;

  /**---------------------------------------------------------------------------------------------------------------------------------/
  /---- RVs produced by this Monte Carlo engine --------------------------------------------------------------------------------------/
  /---- These are free'd when the Monte Carlo engine is destroyed --------------------------------------------------------------------/
  /---------------------------------------------------------------------------------------------------------------------------------**/
  double *variables, *end;

  /**---------------------------------------------------------------------------------------------------------------------------------/
  /---- The density from which RVs are drawn -----------------------------------------------------------------------------------------/
  /---- It does NOT belong to this Monte Carlo engine, but is free'd by the underlying probability space -----------------------------/
  /---------------------------------------------------------------------------------------------------------------------------------**/
  Density* density;

  /**---------------------------------------------------------------------------------------------------------------------------------/
  /---- The class of sampling prototypes, which is used to draw proposals from -------------------------------------------------------/
  /---- It is owned by this Monte Carlo engine and will be free'd, when the Monte Carlo engine is destroyed --------------------------/
  /---- For Metropolis sampling (and up), it must have a method to set transition points ---------------------------------------------/
  /---------------------------------------------------------------------------------------------------------------------------------**/
  Density* sampling_prototype;

  /**---------------------------------------------------------------------------------------------------------------------------------/
  /---- Rejection sampling parameters ------------------------------------------------------------------------------------------------/
  /---------------------------------------------------------------------------------------------------------------------------------**/
  /** The factor, by which the sampling prototype is multiplied. **/
  double rejection_sampling_factor;

} MonteCarloEngine;

typedef struct _ProbabilitySpace{
  /** the dimension of the probability space **/
  int dim;

  /**----------------------------------------------------------------------------------------------------------------------------------/
  /---- RVs controlled by the density                                                                                                --/
  /---- these pointers are owned (and free'd) by this struct                                                                         --/
  /---- the variables, however, are only free'd, if they were automatically created (i.e. they were not registered with this struct) --/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  double **begin, **end;
  
  /** the number of registered RVs **/
  int registered_rvs; 

  /**----------------------------------------------------------------------------------------------------------------------------------/
  /---- A cache for RVs produced by the density ---------------------------------------------------------------------------------------/
  /---- This is only used, if the density can be drawn from directly ------------------------------------------------------------------/
  /---- The cache is free'd, when the probability space is destroyed ------------------------------------------------------------------/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  double *rv_cache_start, *rv_cache_end;

  /**----------------------------------------------------------------------------------------------------------------------------------/
  /---- The density, which governs this probability space -----------------------------------------------------------------------------/
  /---- It needn't be normalized ------------------------------------------------------------------------------------------------------/
  /---- It is destroyed by this probability space, when the probability space is destroyed --------------------------------------------/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  Density *density;

  /**----------------------------------------------------------------------------------------------------------------------------------/
  /---- The Monte Carlo engine, which is used to draw random variables from the density of this probability space ---------------------/
  /---- It is owned (and free'd) by this probability space ----------------------------------------------------------------------------/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  MonteCarloEngine *mce;
} ProbabilitySpace;

typedef struct _GaussianVarianceMatrix{
  /** the dimension of this n*n matrix **/
  int dim;

  /**----------------------------------------------------------------------------------------------------------------------------------/
  /-- the matrix itself ---------------------------------------------------------------------------------------------------------------/
  /-- it is owned by this struct and will be free'd, when the struct is destroyed -----------------------------------------------------/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  double* sigma;

  /**----------------------------------------------------------------------------------------------------------------------------------/
  /-- (LAPACKE) factors used internally to condition the matrix -----------------------------------------------------------------------/
  /-- computed and owned by this struct, it will be free'd when the struct is destroyed -----------------------------------------------/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  double* S;
  /**----------------------------------------------------------------------------------------------------------------------------------/
  /-- (LAPACKE) LU factorization of the matrix sigma ----------------------------------------------------------------------------------/
  /-- it will be free'd, when this struct is destroyed --------------------------------------------------------------------------------/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  double* AF;

  /** determinant of the matrix, calculated from the LU factorization **/
  double det;
} GaussianVarianceMatrix;

typedef struct _GaussianParameters{
  /** dimension of the distribution **/
  int dim;

  /**----------------------------------------------------------------------------------------------------------------------------------/
  /-- expected value of the distribution ----------------------------------------------------------------------------------------------/
  /-- when this struct is destroyed, this array will be free'd if and only if free_mu is true -----------------------------------------/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  const double* mu;
  /** free_mu defaults to true **/
  bool free_mu;

  /**----------------------------------------------------------------------------------------------------------------------------------/
  /-- variance matrix of the distribution ---------------------------------------------------------------------------------------------/
  /-- the flag destroy_gvm indicates whether the matrix should be destroyed and free'd, when this struct is destroyed -----------------/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  GaussianVarianceMatrix* sigma;
  /** destroy_sigma defaults to true **/
  bool destroy_sigma;
} GaussianParameters;

typedef struct _ProductDensityParameters{
  /**----------------------------------------------------------------------------------------------------------------------------------/
  /-- Densities d1 and d2 are currently multiplied and need to have the same dimension ------------------------------------------------/
  /-- CAREFUL: this is actually not a product density in the mathematical sense! ------------------------------------------------------/
  /-- d1 and d2 will be destroyed and free'd, if and only if destroyD1 and destroyD2 are true, respectively. --------------------------/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  Density *d1, *d2;
  /** destroyD1 and destroyD2 default to true **/
  bool destroyD1, destroyD2;
} ProductDensityParameters;

typedef struct _ConditionalDensityParameters{
  /**----------------------------------------------------------------------------------------------------------------------------------/
  /-- The density, which is constrained by conditions ---------------------------------------------------------------------------------/
  /-- It will be destroyed and free'd, if and only if destroyDensity is true ----------------------------------------------------------/
  /----------------------------------------------------------------------------------------------------------------------------------**/
  Density* density;
  /** defaults to true **/
  bool destroyDensity;

  /**-------------------------------------------------------------------------------------------------------------------------------/
  /-- fixed: the fixed values of the conditional density ---------------------------------------------------------------------------/
  /-- parameters: the arguments of the density, which are used as parameters (must be supplied as arguments to evaluate(...) -------/
  /-- positions: the positions of the arguments, which are fixed -------------------------------------------------------------------/
  /---- true: fixed, false: not fixed ----------------------------------------------------------------------------------------------/
  /-- all these arrays will be destroyed, if and only if destroyArrays is true -----------------------------------------------------/
  /-------------------------------------------------------------------------------------------------------------------------------**/
  double* fixed;
  bool* positions;
  bool* parameters;
  /** defaults to true **/
  bool destroyArrays;

  /**-------------------------------------------------------------------------------------------------------------------------------/
  /-- if normalize is true, the density will be normalized first, when it is evaluated ---------------------------------------------/
  /---- this is done by integrating over the density with fixed parameter values ---------------------------------------------------/
  /-------------------------------------------------------------------------------------------------------------------------------**/
  bool normalize;
  
} ConditionalDensityParameters;

typedef struct _ErrorDensityParameters{
  /** a parameter to control the impact of training data **/
  double beta;
  /** the error type **/
  int type;
  /** the dimensions of the inputs, outputs and weights **/
  int dimIN, dimTARGET, dimWEIGHTS;
} ErrorDensityParameters;

typedef struct _ExponentialParameters{
  /** the dimensions of the distribution  **/
  int dim, dimTheta, dimScalarProduct;
  /**----------------------------------------------------------------------------------------------------------------/
  /-- An exponential distribution is a parametrized family of the form ----------------------------------------------/
  /---- d(x|theta) = h(x)*exp(eta(theta)*T(x)-A(theta)) -------------------------------------------------------------/
  /--------------- = h(x)*g(theta)*exp(eta(theta)*T(x)) -------------------------------------------------------------/
  /----------------------------------------------------------------------------------------------------------------**/
  /** when the distribution is destroyed, theta will be free'd, if and only if free_theta is true **/
  double* theta;
  /** defaults to true **/
  bool free_theta;
  /**----------------------------------------------------------------------------------------------------------------/
  /-- these variables are used to store the results of g, A and eta for e specific theta ----------------------------/
  /-- _eta is owned by this struct and free'd, when this distribution is destroyed ----------------------------------/
  /----------------------------------------------------------------------------------------------------------------**/
  double _g, _A; //use nan() to indicate that these are not set
  double *_eta; //in constructor, use 0 to indicate that this is not set, everywhere else use _eta[0] = nan("")
  
  void (*h)(const double* first, const double* last, double* res);
  void (*eta)(const double* first, const double* last, double* res);
  void (*g)(const double* first, const double* last, double* res);
  void (*A)(const double* first, const double* last, double* res);
  void (*T)(const double* first, const double* last, double* res);
} ExponentialParameters;

typedef struct _NetworkStatistics{
  /**-----------------------------------------------------------------------------------------------------------------/
  /-- these arrays belong to this struct -----------------------------------------------------------------------------/
  /-- they will be free'd by destroying it ---------------------------------------------------------------------------/
  /-----------------------------------------------------------------------------------------------------------------**/
  double *mean, *variance;
  double sd, lower_quantile, lower_quartile, median, upper_quartile, upper_quantile;
  bool copy_rvs; // defaults to false
  /**-----------------------------------------------------------------------------------------------------------------/
  /-- a pointer to the output random variables -----------------------------------------------------------------------/
  /-- those are free'd, if and only if copy_rvs is true --------------------------------------------------------------/
  /-----------------------------------------------------------------------------------------------------------------**/
  const double* output_rvs; 
} NetworkStatistics;

#endif
