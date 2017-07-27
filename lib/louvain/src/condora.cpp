// File: condora.cpp
// -- quality functions (for Condorcet-A criterion) source file
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


#include "condora.h"

using namespace std;


CondorA::CondorA(Graph & gr, long double sum):Quality(gr,"A-weighted Condorcet"),sum_se(sum) { 
  n2c.resize(size);

  in.resize(size);

  
  // initialization
  for (int i=0 ; i<size ; i++) {
    n2c[i] = i;
    in[i]  = g.nb_selfloops(i);
  }
}

CondorA::~CondorA() {
  in.clear();
}

long double
CondorA::graph_weighting(Graph *g) {
  long double sum_se = 0.0L;
  
  vector<long double> aux_weights;

  // foreach weight, change Aij to 4Aij/(d(i)+d(i)) - Aii/2d(i) - Ajj/2d(j)
  for (int u=0 ; u < g->nb_nodes ; u++) {
    pair<vector<int>::iterator, vector<long double>::iterator> p = g->neighbors(u);
    int deg = g->nb_neighbors(u);
    for (int i=0 ; i < deg ; i++) {
      int neigh = *(p.first+i);
      long double neigh_w = 0.0L;

      long double aux_neigh_w = 0.0L; // to compute Âij = 2Aij / (d(i)+d(j))
      long double tmp_neigh_w = 0.0L; // to compute (Âii+Âjj)/2 = Aii/2d(i) + Ajj/2d(j)
      
      long double deg_neigh = (long double)(g->nb_neighbors(neigh));
      
      if ((g->weights).size() == 0)
	aux_neigh_w = 2.0L / ((long double)deg + deg_neigh);
      else {
	long double old_neigh = (long double)*(p.second+i);
	aux_neigh_w = 2.0L*old_neigh / ((long double)deg + deg_neigh);
      }

      tmp_neigh_w = (g->nb_selfloops(u)) / (2.0L*(long double)deg) + (g->nb_selfloops(neigh)) / (2.0L*deg_neigh);

      neigh_w = 2.0L*aux_neigh_w - tmp_neigh_w;
      aux_weights.push_back(neigh_w);

      sum_se += tmp_neigh_w - aux_neigh_w;
    }
  }

  g->weights.clear();
  g->weights = aux_weights;

  g->total_weight = 0.0L;

  // Compute total weight
  for (int i=0 ; i < g->nb_nodes ; i++)
    g->total_weight += (long double)(g->weighted_degree(i));

  aux_weights.clear();
  
  return sum_se;
}


long double
CondorA::quality() {
  long double q  = 0.0L;
  long double n  = (long double)g.sum_nodes_w;
  
  for (int i=0 ; i<size ; i++)
    q += in[i];
  
  q += sum_se;
  
  q /= n*n;
  
  return q;
}
