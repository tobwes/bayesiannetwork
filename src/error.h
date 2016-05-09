/* error.h -- This belongs to gneural_network

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

#ifndef ERROR_H
#define ERROR_H

#include "includes.h"
#include "defines.h"
#include "structures.h"
#include "feedforward.h"
#include "network_helpers.h"

extern int NNUM;
extern int NDATA;
extern neuron NEURON[];
extern network NETWORK;

double error(int);

double computeL2SQError(int n); // compute L2 error squared for input X[n] and target Y[n]
double computeLINFError(int n); // compute LINF error for input X[n] and target Y[n]

double specificError(int type, const double* start, const double* target, const double* weights, bool savestate); // error for a specific input/target pair
double computeSpecificError(int type, const double* start, const double* target, const double* weights); // computing the error
double computeSpecificL2SQError(double* output, const double* target); // compute L2 error squared for input/target
double computeSpecificLINFError(double* output, const double* target); // compute LINF error for input/target

#endif
