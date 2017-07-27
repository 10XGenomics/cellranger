// File: condora.h
// -- quality functions (for Condorcet-A criterion) header file
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


#ifndef CONDORA_H
#define CONDORA_H

#include "quality.h"

using namespace std;


class CondorA: public Quality {
 public:
  
  vector<long double> in; // used to compute the quality participation of each community
  long double sum_se; // sum of the nb_selfloops of each vertex

  CondorA(Graph & gr, long double sum);
  ~CondorA();

  // change the weight of each link ij in the graph, from Aij to 4Aij/(d(i)+d(j)) - Aii/2d(i) - Ajj/2d(j)
  // return the result of Sum [Aii/2d(i) + Ajj/2d(j) - 2Aij/(d(i)+d(j))]
  static long double graph_weighting(Graph *g);

  inline void remove(int node, int comm, long double dnodecomm);
  
  inline void insert(int node, int comm, long double dnodecomm);

  inline long double gain(int node, int comm, long double dnodecomm, long double w_degree);

  long double quality();
};


inline void
CondorA::remove(int node, int comm, long double dnodecomm) {
  assert(node>=0 && node<size);

  in[comm] -= 2.0L*dnodecomm + g.nb_selfloops(node);
  
  n2c[node] = -1;
}

inline void
CondorA::insert(int node, int comm, long double dnodecomm) {
  assert(node>=0 && node<size);
  
  in[comm] += 2.0L*dnodecomm + g.nb_selfloops(node);
  
  n2c[node] = comm;
}

inline long double
CondorA::gain(int node, int comm, long double dnc, long double degc) {
  assert(node>=0 && node<size);

  return dnc;
}


#endif // CONDORA_H
