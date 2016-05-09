/* randomize.c -- This belongs to gneural_network

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

// assigns weights randomly for each neuron

#include "randomize.h"

void randomize(void){
  for(int layer=1; layer<NETWORK.num_of_layers; ++layer){
    for(int local_id=0; local_id<NETWORK.num_of_neurons[layer]; ++local_id){
      for(int connection=0; connection<NEURON[NETWORK.neuron_id[layer][local_id]].nw;++connection){
	NEURON[NETWORK.neuron_id[layer][local_id]].w[connection]=WMIN+rnd()*(WMAX-WMIN);
      }
    }
  }
}
