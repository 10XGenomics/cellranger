CC = gcc
CXX = g++
CFLAGS := -Wall -O3 -pipe #-pg
CXXFLAGS := $(CFLAGS)
LOBJS = ssw.o
LCPPOBJS = ssw_cpp.o
PROG = ssw_test
LIB = libssw.so
EXAMPLE = example_c
EXAMPLE_CPP = example_cpp
JAVA_JAR = ssw.jar
JAVA_LIB = libsswjni.so
JAVA_INLCUDES = -I"$(JAVA_HOME)/include" -I"$(JAVA_HOME)/include/linux"
JAVA_OBJ = ssw/Aligner.class ssw/Alignment.class ssw/Example.class

.PHONY: all default java clean

default: $(PROG) $(EXAMPLE) $(EXAMPLE_CPP) $(LIB) 

all: default java

java: $(JAVA_JAR) $(JAVA_LIB)

$(LIB): ssw.c ssw.h
	$(CC) $(CFLAGS) -fPIC -shared -rdynamic -o $@ $<

$(PROG): main.c kseq.h

$(EXAMPLE): example.c

$(PROG) $(EXAMPLE): $(LOBJS)
	$(CC) -o $@ $(filter-out %.h,$^) $(CFLAGS) -lm -lz

$(EXAMPLE_CPP): example.cpp $(LOBJS) $(LCPPOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) -lm -lz

$(JAVA_LIB): sswjni.c ssw.c ssw.h
	$(CC) $(CFLAGS) $(JAVA_INLCUDES) -fPIC -shared -rdynamic -o $@ $< ssw.c 

$(JAVA_JAR): $(JAVA_OBJ)
	jar cvfe $@ ssw.Example $^

%.class: %.java
	javac $<
	
ssw.o: ssw.c ssw.h
	$(CC) -c -o $@ $< $(CFLAGS)

ssw_cpp.o: ssw_cpp.cpp ssw_cpp.h ssw.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	-rm -f $(LOBJS) $(LCPPOBJS) $(PROG) $(LIB) $(EXAMPLE) $(EXAMPLE_CPP) $(JAVA_LIB) $(JAVA_JAR) $(JAVA_OBJ) *~ 
