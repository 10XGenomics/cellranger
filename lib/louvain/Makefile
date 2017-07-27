#!/bin/bash

CXX=g++
CXXFLAGS= -ansi -O5 -Wall
DIRSRC= ./src/
EXEC=louvain convert hierarchy matrix
OBJ1= $(DIRSRC)graph_binary.o $(DIRSRC)louvain.o $(DIRSRC)quality.o $(DIRSRC)modularity.o $(DIRSRC)zahn.o $(DIRSRC)owzad.o $(DIRSRC)goldberg.o $(DIRSRC)condora.o $(DIRSRC)devind.o $(DIRSRC)devuni.o $(DIRSRC)dp.o $(DIRSRC)shimalik.o $(DIRSRC)balmod.o
OBJ2= $(DIRSRC)graph.o 

all: $(EXEC)

louvain : $(OBJ1) $(DIRSRC)main_louvain.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

convert : $(OBJ2) $(DIRSRC)main_convert.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

hierarchy : $(DIRSRC)main_hierarchy.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

matrix : $(DIRSRC)main_matrix.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CXX) -o  $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(DIRSRC)*.o

mrproper: clean
	rm -f *~ $(EXEC)
