// File: shimalik.cpp
// -- quality functions (for Shi-Malik criterion) source file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// And based on the article "A Generalized and Adaptive Method for Community Detection"
// Copyright (C) 2014 R. Campigotto, P. Conde Céspedes, J.-L. Guillaume
//
// This file is part of Louvain algorithm.
// 
// Louvain algorithm is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Louvain algorithm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : July 2014
//-----------------------------------------------------------------------------
// see README.txt for more details


#include "shimalik.h"

using namespace std;


ShiMalik::ShiMalik(Graph & gr, int kappa_min):Quality(gr,"Shi-Malik (with kmin=" + to_string(kappa_min) + ")"), kappa(size),kmin(kappa_min) {
  n2c.resize(size);

  in.resize(size);
  tot.resize(size);

  // initialization
  for (int i=0 ; i<size ; i++) {
    n2c[i] = i;
    in[i]  = g.nb_selfloops(i);
    tot[i] = g.weighted_degree(i);
  }
}

ShiMalik::~ShiMalik() {
  in.clear();
  tot.clear();
}

long double
ShiMalik::quality() {
  long double q  = 0.0L;
  long double n = (long double)g.sum_nodes_w;

  for (int i=0 ; i<size ; i++) {
    if (tot[i] > 0.0L)
      q += in[i] / tot[i];
  }
  
  q -= (long double)kappa;

  q /= n;
  
  return q;
}
