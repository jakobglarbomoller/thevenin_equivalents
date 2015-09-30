# makefile for gtc
#

CXX?=g++

LOCAL_LIBS= -lamd -lcxsparse -lboost_system -lboost_timer -lglog -lgflags
LOCAL_OPTS= -O3 -std=gnu++11 -I/usr/include/suitesparse -I/usr/include/boost -fopenmp


EXECS=main

WDIR := ${CURDIR}

.PHONY: all
all: main


compute_gtc.o: compute_gtc.cpp compute_gtc.hpp
	${CXX} -c $< -o $@ ${LOCAL_OPTS}


main: main.cpp compute_gtc.o
	${CXX} $< -o $@ ${LOCAL_OPTS} $(filter %.o,$^) ${LOCAL_LIBS}


.PHONY: clean
clean:
	rm -f ${EXECS} *.o core* *~ 
