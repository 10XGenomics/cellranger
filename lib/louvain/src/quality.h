// File: quality.h
// -- quality functions header file
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


#ifndef QUALITY_H
#define QUALITY_H

#include <sstream>

#include "graph_binary.h"

using namespace std;


class Quality {
 public:
  
  Graph & g; // network to compute communities for
  int size; // nummber of nodes in the network and size of all vectors
  string name;
  
  vector<int> n2c; // community to which each node belongs
 Quality(Graph &gr, const std::string& n):g(gr),size(g.nb_nodes),name(n){}
  
  virtual ~Quality();
  
  // remove the node from its current community with which it has dnodecomm links
  virtual void remove(int node, int comm, long double dnodecomm)=0;
  
  // insert the node in comm with which it shares dnodecomm links
  virtual void insert(int node, int comm, long double dnodecomm)=0;
  
  // compute the gain of quality by adding node to comm
  virtual long double gain(int node, int comm, long double dnodecomm, long double w_degree)=0;
  
  // compute the quality of the current partition
  virtual long double quality()=0;
};

template<class T>
std::string to_string(T number){
  std::ostringstream oss;
  oss << number;
  return oss.str();
}

#endif // QUALITY_H
