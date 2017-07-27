// File: louvain.cpp
// -- community detection source file
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
// see readme.txt for more details

#include "louvain.h"

using namespace std;


Louvain::Louvain(int nbp, long double epsq, Quality* q) {
  qual = q;

  neigh_weight.resize(qual->size,-1);
  neigh_pos.resize(qual->size);
  neigh_last = 0;

  nb_pass = nbp;
  eps_impr = epsq;
}

void
Louvain::init_partition(char * filename) {
  ifstream finput;
  finput.open(filename,fstream::in);

  // read partition
  while (!finput.eof()) {
    int node, comm;
    finput >> node >> comm;
    
    if (finput) {
      int old_comm = qual->n2c[node];
      neigh_comm(node);

      qual->remove(node, old_comm, neigh_weight[old_comm]);
      
      int i=0;
      for (i=0 ; i<neigh_last ; i++) {
	int best_comm = neigh_pos[i];
	long double best_nblinks = neigh_weight[neigh_pos[i]];
	if (best_comm==comm) {
	  qual->insert(node, best_comm, best_nblinks);
	  break;
	}
      }
      if (i==neigh_last)
	qual->insert(node, comm, 0);
    }
  }
  finput.close();
}

void
Louvain::neigh_comm(int node) {
  for (int i=0 ; i<neigh_last ; i++)
    neigh_weight[neigh_pos[i]]=-1;
  
  neigh_last = 0;

  pair<vector<int>::iterator, vector<long double>::iterator> p = (qual->g).neighbors(node);
  int deg = (qual->g).nb_neighbors(node);

  neigh_pos[0] = qual->n2c[node];
  neigh_weight[neigh_pos[0]] = 0;
  neigh_last = 1;

  for (int i=0 ; i<deg ; i++) {
    int neigh  = *(p.first+i);
    int neigh_comm = qual->n2c[neigh];
    long double neigh_w = ((qual->g).weights.size()==0)?1.0L:*(p.second+i);
    
    if (neigh!=node) {
      if (neigh_weight[neigh_comm]==-1) {
	neigh_weight[neigh_comm] = 0.0L;
	neigh_pos[neigh_last++] = neigh_comm;
      }
      neigh_weight[neigh_comm] += neigh_w;
    }
  }
}

void
Louvain::partition2graph() {
  vector<int> renumber(qual->size, -1);
  for (int node=0 ; node<qual->size ; node++) {
    renumber[qual->n2c[node]]++;
  }

  int end=0;
  for (int i=0 ; i< qual->size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=end++;

  for (int i=0 ; i< qual->size ; i++) {
    pair<vector<int>::iterator, vector<long double>::iterator> p = (qual->g).neighbors(i);

    int deg = (qual->g).nb_neighbors(i);
    for (int j=0 ; j<deg ; j++) {
      int neigh = *(p.first+j);
      cout << renumber[qual->n2c[i]] << " " << renumber[qual->n2c[neigh]] << endl;
    }
  }
}

void
Louvain::display_partition() {
  vector<int> renumber(qual->size, -1);
  for (int node=0 ; node < qual->size ; node++) {
    renumber[qual->n2c[node]]++;
  }

  int end=0;
  for (int i=0 ; i < qual->size ; i++)
    if (renumber[i]!=-1)
      renumber[i] = end++;

  for (int i=0 ; i < qual->size ; i++)
    cout << i << " " << renumber[qual->n2c[i]] << endl;
}

Graph
Louvain::partition2graph_binary() {
  // Renumber communities
  vector<int> renumber(qual->size, -1);
  for (int node=0 ; node < qual->size ; node++)
    renumber[qual->n2c[node]]++;

  int last=0;
  for (int i=0 ; i < qual->size ; i++) {
    if (renumber[i]!=-1)
      renumber[i] = last++;
  }
  
  // Compute communities
  vector<vector<int> > comm_nodes(last);
  vector<int> comm_weight(last, 0);
  
  for (int node = 0 ; node < (qual->size) ; node++) {
    comm_nodes[renumber[qual->n2c[node]]].push_back(node);
    comm_weight[renumber[qual->n2c[node]]] += (qual->g).nodes_w[node];
  }

  // Compute weighted graph
  Graph g2;
  int nbc = comm_nodes.size();

  g2.nb_nodes = comm_nodes.size();
  g2.degrees.resize(nbc);
  g2.nodes_w.resize(nbc);
  
  for (int comm=0 ; comm<nbc ; comm++) {
    map<int,long double> m;
    map<int,long double>::iterator it;

    int size_c = comm_nodes[comm].size();

    g2.assign_weight(comm, comm_weight[comm]);

    for (int node=0 ; node<size_c ; node++) {
      pair<vector<int>::iterator, vector<long double>::iterator> p = (qual->g).neighbors(comm_nodes[comm][node]);
      int deg = (qual->g).nb_neighbors(comm_nodes[comm][node]);
      for (int i=0 ; i<deg ; i++) {
	int neigh = *(p.first+i);
	int neigh_comm = renumber[qual->n2c[neigh]];
	long double neigh_weight = ((qual->g).weights.size()==0)?1.0L:*(p.second+i);

	it = m.find(neigh_comm);
	if (it==m.end())
	  m.insert(make_pair(neigh_comm, neigh_weight));
	else
	  it->second += neigh_weight;
      }
    }

    g2.degrees[comm] = (comm==0)?m.size():g2.degrees[comm-1]+m.size();
    g2.nb_links += m.size();

    for (it = m.begin() ; it!=m.end() ; it++) {
      g2.total_weight += it->second;
      g2.links.push_back(it->first);
      g2.weights.push_back(it->second);
    }
  }

  return g2;
}

bool
Louvain::one_level() {
  bool improvement=false ;
  int nb_moves;
  int nb_pass_done = 0;
  long double new_qual = qual->quality();
  long double cur_qual = new_qual;

  vector<int> random_order(qual->size);
  for (int i=0 ; i < qual->size ; i++)
    random_order[i]=i;
  for (int i=0 ; i < qual->size-1 ; i++) {
    int rand_pos = rand()%(qual->size-i)+i;
    int tmp = random_order[i];
    random_order[i] = random_order[rand_pos];
    random_order[rand_pos] = tmp;
  }

  // repeat while 
  //   there is an improvement of quality
  //   or there is an improvement of quality greater than a given epsilon 
  //   or a predefined number of pass have been done
  do {
    cur_qual = new_qual;
    nb_moves = 0;
    nb_pass_done++;

    // for each node: remove the node from its community and insert it in the best community
    for (int node_tmp = 0 ; node_tmp < qual->size ; node_tmp++) {
      int node = random_order[node_tmp];
      int node_comm = qual->n2c[node];
      long double w_degree = (qual->g).weighted_degree(node);

      // computation of all neighboring communities of current node
      neigh_comm(node);
      // remove node from its current community
      qual->remove(node, node_comm, neigh_weight[node_comm]);

      // compute the nearest community for node
      // default choice for future insertion is the former community
      int best_comm = node_comm;
      long double best_nblinks  = 0.0L;
      long double best_increase = 0.0L;
      for (int i=0 ; i<neigh_last ; i++) {
	long double increase = qual->gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree);
	if (increase>best_increase) {
	  best_comm = neigh_pos[i];
	  best_nblinks = neigh_weight[neigh_pos[i]];
	  best_increase = increase;
	}
      }

      // insert node in the nearest community
      qual->insert(node, best_comm, best_nblinks);
     
      if (best_comm!=node_comm)
	nb_moves++;
    }

    new_qual = qual->quality();
    
    if (nb_moves>0)
      improvement=true;

  } while (nb_moves>0 && new_qual-cur_qual > eps_impr);

  return improvement;
}
