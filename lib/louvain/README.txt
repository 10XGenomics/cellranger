-----------------------------------------------------------------------------

Multi-criteria community detection
Version 0.3 - not compatible with the previous version (see below).

Based on the article "Fast unfolding of community hierarchies in large networks"
Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre

And based on the article "A Generalized and Adaptive Method for Community Detection"
Copyright (C) 2014 R. Campigotto, P. Conde Céspedes, J.-L. Guillaume

This file is part of Louvain algorithm.

Louvain algorithm is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Louvain algorithm is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.

-----------------------------------------------------------------------------

Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto
Email    : jean-loup.guillaume@lip6.fr
Location : Paris, France
Time	 : July 2014

-----------------------------------------------------------------------------

Disclaimer:
If you find a bug, please send a bug report to jean-loup.guillaume@lip6.fr
including if necessary the input file and the parameters that caused the bug.
You can also send me any comment or suggestion about the program.

Note that the program is expecting a friendly use and therefore does not make
much verifications about the arguments.

-----------------------------------------------------------------------------


This package offers a set of functions to use in order to compute 
communities on graphs weighted or unweighted. A typical sequence of 
actions is:

1. Conversion from a text format (each line contains a couple "src dest")
./convert -i graph.txt -o graph.bin
This program can also be used to convert weighted graphs (each line contain
a triple "src dest w") using -w option:
./convert -i graph.txt -o graph.bin -w graph.weights
Finally, nodes can be renumbered from 0 to nb_nodes-1 using -r option
(less space wasted in some cases):
./convert -i graph.txt -o graph.bin -r labelings_connection_file.txt


2. Computes communities with a specified quality function and displays hierarchical tree:
./louvain graph.bin -l -1 -v -q id_qual > graph.tree

To ensure a faster computation (with a loss of quality), one can use
the -e option to specify that the program must stop if the increase of
modularity is below epsilon for a given iteration or pass:
./louvain graph.bin -l -1 -q id_qual -e 0.001 > graph.tree

The program can deal with weighted networks using -w option:
./louvain graph.bin -l -1 -q id_qual -w graph.weights > graph.tree
In this specific case, the convertion step must also use the -w option.

The program can also start with any given partition using -p option
./louvain graph.bin -q id_qual -p graph.part -v


3. Displays information on the tree structure (number of hierarchical
levels and nodes per level):
./hierarchy graph.tree

Displays the belonging of nodes to communities for a given level of the tree:
./hierarchy graph.tree -l 2 > graph_node2comm_level2

Displays the belonging of nodes to communities for the last level of the tree:
./hierarchy graph.tree -m > graph_node2comm_lastlevel


4. To display the X relational matrix for a given level of the tree:
./matrix graph.tree -l 2 > graph_X_level2

Displays the X relational matrix for the last level of the tree:
./matrix graph.tree -m > graph_node2comm_lastlevel  



-----------------------------------------------------------------------------

Known bugs or restrictions:
- the number of nodes is stored on int (4 bytes) and the number of links on unsigned long long (8 bytes).
- the program is not suited for 32 bits Windows systems.

-----------------------------------------------------------------------------

Version history:
The following modifications have been made from version 0.1 and 0.2:
- weights are now stored using long double (integer in V0.1 and float in V0.2)
- degrees are stored on 8 bytes (unsigned long long) allowing large graphs to be decomposed
- weights are stored in a separate file, which allows disk usage reduction if
  different weights are to be used on the same topology
- any given partition can be used as a seed for the algorithm rather than just
  the trivial partition where each node belongs to its own community
- initial network can contain loops if network is considered weighted
- graph is not renumbered by default in the convert program
- an optional verbose mode has been added and the program is silent by default
- some portions of the code have been c++ improved (type * -> vector<type>)
These modifications imply that any binary graph file created with the previous 
versions of the code is not comptabile with this version. You must therefore regenerate all the binary files.

Version 0.1 and 0.2:
- initial community detection algorithm

