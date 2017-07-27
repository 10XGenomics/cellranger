// File: main_louvain.cpp
// -- community detection, sample main file
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


#include "graph_binary.h"
#include "louvain.h"
#include <unistd.h>

#include "modularity.h"
#include "zahn.h"
#include "owzad.h"
#include "goldberg.h"
#include "condora.h"
#include "devind.h"
#include "devuni.h"
#include "dp.h"
#include "shimalik.h"
#include "balmod.h"


using namespace std;


char *filename = NULL;
char *filename_w = NULL;
char *filename_part = NULL;
int type = UNWEIGHTED;

int nb_pass = 0;
long double precision = 0.000001L;
int display_level = -2;

unsigned short id_qual = 0;

long double alpha = 0.5L;
int kmin = 1;

long double sum_se = 0.0L;
long double sum_sq = 0.0L;

long double max_w = 1.0L;

Quality *q;

bool verbose = false;

void
usage(char *prog_name, const char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " input_file [-q id_qual] [-c alpha] [-k min] [-w weight_file] [-p part_file] [-e epsilon] [-l display_level] [-v] [-h]" << endl << endl;
  cerr << "input_file: file containing the graph to decompose in communities" << endl;

  cerr << "-q id\tthe quality function used to compute partition of the graph (modularity is chosen by default):" << endl << endl;

  cerr << "\tid = 0\t -> the classical Newman-Girvan criterion (also called \"Modularity\")" << endl;
  cerr << "\tid = 1\t -> the Zahn-Condorcet criterion" << endl;
  cerr << "\tid = 2\t -> the Owsinski-Zadrozny criterion (you should specify the value of the parameter with option -c)" << endl;
  cerr << "\tid = 3\t -> the Goldberg Density criterion" << endl;
  cerr << "\tid = 4\t -> the A-weighted Condorcet criterion" << endl;
  cerr << "\tid = 5\t -> the Deviation to Indetermination criterion" << endl;
  cerr << "\tid = 6\t -> the Deviation to Uniformity criterion" << endl;
  cerr << "\tid = 7\t -> the Profile Difference criterion" << endl;
  cerr << "\tid = 8\t -> the Shi-Malik criterion (you should specify the value of kappa_min with option -k)" << endl;
  cerr << "\tid = 9\t -> the Balanced Modularity criterion" << endl;

  cerr << endl;

  cerr << "-c al\tthe parameter for the Owsinski-Zadrozny quality function (between 0.0 and 1.0: 0.5 is chosen by default)" << endl;
  cerr << "-k min\tthe kappa_min value (for Shi-Malik quality function) (it must be > 0: 1 is chosen by default)" << endl;

  cerr << endl;

  cerr << "-w file\tread the graph as a weighted one (weights are set to 1 otherwise)" << endl;
  cerr << "-p file\tstart the computation with a given partition instead of the trivial partition" << endl;
  cerr << "\tfile must contain lines \"node community\"" << endl;
  cerr << "-e eps\ta given pass stops when the quality is increased by less than epsilon" << endl;
  cerr << "-l k\tdisplays the graph of level k rather than the hierachical structure" << endl;
  cerr << "\tif k=-1 then displays the hierarchical structure rather than the graph at a given level" << endl;
  cerr << "-v\tverbose mode: gives computation time, information about the hierarchy and quality" << endl;
  cerr << "-h\tshow this usage message" << endl;

  exit(0);
}

void
parse_args(int argc, char **argv) {
  if (argc<2)
    usage(argv[0], "Bad arguments number\n");

  for (int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
      case 'w':
	type = WEIGHTED;
        filename_w = argv[i+1];
	i++;
	break;
      case 'q':
	id_qual = (unsigned short)atoi(argv[i+1]);
	i++;
	break;
      case 'c':
	alpha = atof(argv[i+1]);
	i++;
	break;
      case 'k':
	kmin = atoi(argv[i+1]);
	i++;
	break;
      case 'p':
        filename_part = argv[i+1];
	i++;
	break;
      case 'e':
	precision = atof(argv[i+1]);
	i++;
	break;
      case 'l':
	display_level = atoi(argv[i+1]);
	i++;
	break;
      case 'v':
	verbose = true;
	break;
      case 'h':
	usage(argv[0], "");
	break;
      default:
	usage(argv[0], "Unknown option\n");
      }
    } else {
      if (filename==NULL)
        filename = argv[i];
      else
        usage(argv[0], "More than one filename\n");
    }
  }
  if (filename == NULL)
    usage(argv[0], "No input file has been provided\n");
}

void
display_time(const char *str) {
  time_t rawtime;
  time ( &rawtime );
  cerr << str << ": " << ctime (&rawtime);
}

void
init_quality(Graph *g, unsigned short nbc) {

  if (nbc > 0)
    delete q;

  switch (id_qual) {
  case 0:
    q = new Modularity(*g);
    break;
  case 1:
    if (nbc == 0)
      max_w = g->max_weight();
    q = new Zahn(*g, max_w);
    break;
  case 2:
    if (nbc == 0)
      max_w = g->max_weight();
    if (alpha <= 0. || alpha >= 1.0)
      alpha = 0.5;
    q = new OwZad(*g, alpha, max_w);
    break;
  case 3:
    if (nbc == 0)
      max_w = g->max_weight();
    q = new Goldberg(*g, max_w);
    break;
  case 4:
    if (nbc == 0) {
      g->add_selfloops();
      sum_se = CondorA::graph_weighting(g);
    }
    q = new CondorA(*g, sum_se);
    break;
  case 5:
    q = new DevInd(*g);
    break;
  case 6:
    q = new DevUni(*g);
    break;
  case 7:
    if (nbc == 0) {
      max_w = g->max_weight();
      sum_sq = DP::graph_weighting(g);
    }
    q = new DP(*g, sum_sq, max_w);
    break;
  case 8:
    if (kmin < 1)
      kmin = 1;
    q = new ShiMalik(*g, kmin);
    break;
  case 9:
    if (nbc == 0)
      max_w = g->max_weight();
    q = new BalMod(*g, max_w);
    break;
  default:
    q = new Modularity(*g);
    break;
  }
}


int
main(int argc, char **argv) {

  srand(0x00C0FFEE);
  
  parse_args(argc, argv);
  
  time_t time_begin, time_end;
  time(&time_begin);
  
  unsigned short nb_calls = 0;
  
  if (verbose) 
    display_time("Begin");
  
  Graph g(filename, filename_w, type);
  init_quality(&g, nb_calls);
  nb_calls++;
  
  if (verbose)
    cerr << endl << "Computation of communities with the " << q->name << " quality function" << endl << endl;
  
  Louvain c(-1, precision, q);
  if (filename_part!=NULL)
    c.init_partition(filename_part);
  
  bool improvement = true;
  
  long double quality = (c.qual)->quality();
  long double new_qual;
  
  int level = 0;
  
  do {
    if (verbose) {
      cerr << "level " << level << ":\n";
      display_time("  start computation");
      cerr << "  network size: " 
	   << (c.qual)->g.nb_nodes << " nodes, " 
	   << (c.qual)->g.nb_links << " links, "
	   << (c.qual)->g.total_weight << " weight" << endl;
    }
    
    improvement = c.one_level();
    new_qual = (c.qual)->quality();
    
    if (++level==display_level)
      (c.qual)->g.display();
    if (display_level==-1)
      c.display_partition();
    
    g = c.partition2graph_binary();
    init_quality(&g, nb_calls);
    nb_calls++;
    
    c = Louvain(-1, precision, q);
    
    if (verbose)
      cerr << "  quality increased from " << quality << " to " << new_qual << endl;
    
    //quality = (c.qual)->quality();
    quality = new_qual;
    
    if (verbose)
      display_time("  end computation");
    
    if (filename_part!=NULL && level==1) // do at least one more computation if partition is provided
      improvement=true;
  } while(improvement);
  
  time(&time_end);
  if (verbose) {
    display_time("End");
    cerr << "Total duration: " << (time_end-time_begin) << " sec" << endl;
  }
  cerr << new_qual << endl;
  
  delete q;
}
