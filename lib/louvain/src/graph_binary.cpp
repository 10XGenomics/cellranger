// File: graph_binary.cpp
// -- graph handling source
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


#include <fstream>
#include "graph_binary.h"


Graph::Graph() {
  nb_nodes = 0;
  nb_links = 0ULL;

  total_weight = 0.0L;
  sum_nodes_w = 0;
}

Graph::Graph(char *filename, char *filename_w, int type) {
  ifstream finput;
  finput.open(filename,fstream::in | fstream::binary);
  if (finput.is_open() != true) {
    cerr << "The file " << filename << " does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  // Read number of nodes on 4 bytes
  finput.read((char *)&nb_nodes, sizeof(int));
  if (finput.rdstate() != ios::goodbit) {
    cerr << "The file " << filename << " is not a valid graph" << endl;
    exit(EXIT_FAILURE);
  }
  
  // Read cumulative degree sequence: 8 bytes for each node
  // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
  degrees.resize(nb_nodes);
  finput.read((char *)&degrees[0], nb_nodes*sizeof(unsigned long long));

  // Read links: 4 bytes for each link (each link is counted twice)
  nb_links = degrees[nb_nodes-1];
  links.resize(nb_links);
  finput.read((char *)(&links[0]), nb_links*sizeof(int));

  // IF WEIGHTED, read weights: 10 bytes for each link (each link is counted twice)
  weights.resize(0);
  total_weight = 0.0L;
  if (type==WEIGHTED) {
    ifstream finput_w;
    finput_w.open(filename_w,fstream::in | fstream::binary);
    if (finput_w.is_open() != true) {
      cerr << "The file " << filename_w << " does not exist" << filename << endl;
      exit(EXIT_FAILURE);
    }

    weights.resize(nb_links);
    finput_w.read((char *)(&weights[0]), nb_links*sizeof(long double));
    if (finput_w.rdstate() != ios::goodbit) {
      cerr << "The file " << filename_w << " does not correspond to valid weights for the graph" << filename << endl;
      exit(EXIT_FAILURE);
    }
  }

  // Compute total weight
  for (int i=0 ; i<nb_nodes ; i++)
    total_weight += (long double)weighted_degree(i);

  nodes_w.assign(nb_nodes, 1);
  sum_nodes_w = nb_nodes;
}

long double
Graph::max_weight() {
  long double max = 1.0L;

  if (weights.size()!=0)
    max = *max_element(weights.begin(),weights.end());
  
  return max;
}

void
Graph::assign_weight(int node, int weight) {
  sum_nodes_w -= nodes_w[node];

  nodes_w[node] = weight;

  sum_nodes_w += weight;
}

void
Graph::add_selfloops() {
  vector<unsigned long long> aux_deg;
  vector<int> aux_links;

  unsigned long long sum_d = 0ULL;

  for (int u=0 ; u < nb_nodes ; u++) {
    pair<vector<int>::iterator, vector<long double>::iterator> p = neighbors(u);
    int deg = nb_neighbors(u);

    for (int i=0 ; i < deg ; i++) {
      int neigh = *(p.first+i);
      aux_links.push_back(neigh);
    }

    sum_d += (unsigned long long)deg;

    if (nb_selfloops(u) == 0.0L) {
      aux_links.push_back(u); // add a selfloop
      sum_d += 1ULL;
    }

    aux_deg.push_back(sum_d); // add the (new) degree of vertex u
  }

  links = aux_links;
  degrees = aux_deg;
  
  nb_links += (unsigned long long)nb_nodes;
}

void
Graph::display() {
  for (int node=0 ; node<nb_nodes ; node++) {
    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
    cout << node << ":" ;
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      if (true) {
	if (weights.size()!=0)
	  cout << " (" << *(p.first+i) << " " << *(p.second+i) << ")";
	else
	  cout << " " << *(p.first+i);
      }
    }
    cout << endl;
  }
}

void
Graph::display_reverse() {
  for (int node=0 ; node<nb_nodes ; node++) {
    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      if (node>*(p.first+i)) {
	if (weights.size()!=0)
	  cout << *(p.first+i) << " " << node << " " << *(p.second+i) << endl;
	else
	  cout << *(p.first+i) << " " << node << endl;
      }
    }
  }
}

bool
Graph::check_symmetry() {
  int error = 0;
  for (int node=0 ; node<nb_nodes ; node++) {
    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      int neigh = *(p.first+i);
      long double weight = *(p.second+i);

      pair<vector<int>::iterator, vector<long double>::iterator > p_neigh = neighbors(neigh);
      for (int j=0 ; j<nb_neighbors(neigh) ; j++) {
	int neigh_neigh = *(p_neigh.first+j);
	long double neigh_weight = *(p_neigh.second+j);

	if (node==neigh_neigh && weight!=neigh_weight) {
	  cout << node << " " << neigh << " " << weight << " " << neigh_weight << endl;
	  if (error++==10)
	    exit(0);
	}
      }
    }
  }
  return (error==0);
}

void
Graph::display_binary(char *outfile) {
  ofstream foutput;
  foutput.open(outfile ,fstream::out | fstream::binary);

  foutput.write((char *)(&nb_nodes),sizeof(int));
  foutput.write((char *)(&degrees[0]),sizeof(unsigned long long)*nb_nodes);
  foutput.write((char *)(&links[0]),sizeof(int)*nb_links);
}

