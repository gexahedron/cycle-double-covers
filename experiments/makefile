SHELL = /bin/sh
#CC = clang
CC = g++
CFLAGS = -O3 --std=c++1z -fstack-protector-all
INCLUDES += -I .

all :
	${CC} $(CFLAGS) $(INCLUDES) \
	graph.cpp args.cpp \
	util/flows.cpp \
	experiments/cycles.cpp \
	experiments/preimages.cpp \
	experiments/o5cdc.cpp \
	experiments/o6c4c.cpp \
	experiments/unit_vector_flows.cpp \
	experiments/petersen_colouring.cpp \
	main.cpp -o run_experiments
