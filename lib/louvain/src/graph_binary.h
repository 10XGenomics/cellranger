// File: graph_binary.h
// -- graph handling header file
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


#ifndef GRAPH_H
#define GRAPH_H

#include <assert.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;


class Graph {
 public:
  int nb_nodes;
  unsigned long long nb_links;

  long double total_weight;
  int sum_nodes_w;

  vector<unsigned long long> degrees;
  vector<int> links;
  vector<long double> weights;

  vector<int> nodes_w;

  Graph();
  
  // binary file format is
  // 4 bytes for the number of nodes in the graph
  // 8*(nb_nodes) bytes for the cumulative degree for each node:
  //    deg(0)=degrees[0]
  //    deg(k)=degrees[k]-degrees[k-1]
  // 4*(sum_degrees) bytes for the links
  // IF WEIGHTED, 10*(sum_degrees) bytes for the weights in a separate file
  Graph(char *filename, char *filename_w, int type);

  // return the biggest weight of links in the graph
  long double max_weight();
  
  // assign a weight to a node (needed after the first level)
  void assign_weight(int node, int weight);

  // add selfloop to each vertex in the graph
  void add_selfloops();

  void display(void);
  void display_reverse(void);
  void display_binary(char *outfile);
  bool check_symmetry();


  // return the number of neighbors (degree) of the node
  inline int nb_neighbors(int node);

  // return the number of self loops of the node
  inline long double nb_selfloops(int node);

  // return the weighted degree of the node
  inline long double weighted_degree(int node);

  // return pointers to the first neighbor and first weight of the node
  inline pair<vector<int>::iterator, vector<long double>::iterator > neighbors(int node);
};


inline int
Graph::nb_neighbors(int node) {
  assert(node>=0 && node<nb_nodes);

  if (node==0)
    return degrees[0];
  else
    return (int)(degrees[node]-degrees[node-1]);
}

inline long double
Graph::nb_selfloops(int node) {
  assert(node>=0 && node<nb_nodes);

  pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
  for (int i=0 ; i<nb_neighbors(node) ; i++) {
    if (*(p.first+i)==node) {
      if (weights.size()!=0)
	return (long double)*(p.second+i);
      else 
	return 1.0L;
    }
  }
  return 0.0L;
}

inline long double
Graph::weighted_degree(int node) {
  assert(node>=0 && node<nb_nodes);
  
  if (weights.size()==0)
    return (long double)nb_neighbors(node);
  else {
    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
    long double res = 0.0L;
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      res += (long double)*(p.second+i);
    }
    return res;
  }
}

inline pair<vector<int>::iterator, vector<long double>::iterator >
Graph::neighbors(int node) {
  assert(node>=0 && node<nb_nodes);
  
  if (node==0)
    return make_pair(links.begin(), weights.begin());
  else if (weights.size()!=0)
    return make_pair(links.begin()+degrees[node-1], weights.begin()+degrees[node-1]);
  else
    return make_pair(links.begin()+degrees[node-1], weights.begin());
}


#endif // GRAPH_H
