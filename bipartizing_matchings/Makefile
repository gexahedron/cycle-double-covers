CXX ?= g++
CXXFLAGS += -std=c++11 -O3
INCLUDES += -I ../util
UTIL += ../util/util.cpp
PROGRAM_NAME = bipartizing_matchings

all :
	${CXX} $(CXXFLAGS) $(INCLUDES) $(PROGRAM_NAME).cpp $(UTIL) -o $(PROGRAM_NAME)
