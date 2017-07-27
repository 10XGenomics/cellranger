// File: modularity.h
// -- quality functions (for Modularity criterion) header file
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


#ifndef MODULARITY_H
#define MODULARITY_H

#include "quality.h"

using namespace std;


class Modularity: public Quality {
 public:

  vector<long double> in, tot; // used to compute the quality participation of each community

  Modularity(Graph & gr);
  ~Modularity();

  inline void remove(int node, int comm, long double dnodecomm);

  inline void insert(int node, int comm, long double dnodecomm);

  inline long double gain(int node, int comm, long double dnodecomm, long double w_degree);

  long double quality();
};


inline void
Modularity::remove(int node, int comm, long double dnodecomm) {
  assert(node>=0 && node<size);

  in[comm]  -= 2.0L*dnodecomm + g.nb_selfloops(node);
  tot[comm] -= g.weighted_degree(node);
  
  n2c[node] = -1;
}

inline void
Modularity::insert(int node, int comm, long double dnodecomm) {
  assert(node>=0 && node<size);
  
  in[comm]  += 2.0L*dnodecomm + g.nb_selfloops(node);
  tot[comm] += g.weighted_degree(node);
  
  n2c[node] = comm;
}

inline long double
Modularity::gain(int node, int comm, long double dnc, long double degc) {
  assert(node>=0 && node<size);
  
  long double totc = tot[comm];
  long double m2   = g.total_weight;
  
  return (dnc - totc*degc/m2);
}


#endif // MODULARITY_H
