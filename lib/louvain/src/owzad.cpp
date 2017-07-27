// File: owzad.cpp
// -- quality functions (for Owsinski-Zadrozny criterion) source file
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


#include "owzad.h"

using namespace std;


OwZad::OwZad(Graph & gr, long double al, long double max_w):Quality(gr,"Owsinski-Zadrozny (with alpha=" + to_string(al) + ")"), alpha(al), max(max_w) {  
  n2c.resize(size);

  in.resize(size);
  w.resize(size);

  
  // initialization
  for (int i=0 ; i<size ; i++) {
    n2c[i] = i;
    in[i]  = g.nb_selfloops(i);
    w[i]   = g.nodes_w[i];
  }

}

OwZad::~OwZad() {
  in.clear();
  w.clear();
}

long double
OwZad::quality() {
  long double q  = 0.0L;
  long double m2 = g.total_weight;
  long double n  = (long double)g.sum_nodes_w;

  for (int i=0 ; i < size ; i++) {
    long double wc = (long double)w[i];
    if (wc > 0.0L)
      q += in[i] - alpha*(wc*wc*max);
  }
  
  q += alpha * (n*n*max - m2);
  
  q /= n*n*max;
  
  return q;
}
