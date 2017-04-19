# Feb 2015, Stephen Becker
# Run this, or run the Matlab script
# (on Windows, the Matlab script is probably easier)
# Assumes you have gcc


MATLAB = $(shell matlab -e | sed -n 's/MATLAB=//p')
CC = $(MATLAB)/bin/mex
#CC = mex  # Can be mistaken for pdftex
FLAGS = -largeArrayDims -lm
FLAGS += -O -g
#FLAGS += "CFLAGS='$$CFLAGS -Wall'"
#-Wall -Wno-uninitialized
# If you want, define or undefine the debug flag
#FLAGS += -DDEBUG
FLAGS += -UDEBUG

SRC_DIR=../src

INCLUDES = -I$(SRC_DIR)

LBFGSB  = $(SRC_DIR)/lbfgsb.c $(SRC_DIR)/linesearch.c \
		  $(SRC_DIR)/subalgorithms.c $(SRC_DIR)/print.c

LINPACK = $(SRC_DIR)/linpack.c

BLAS 	= $(SRC_DIR)/miniCBLAS.c
#CFLAGS += -D_USE_OPTIMIZED_BLAS -lblas
#CFLAGS += -D_USE_OPTIMIZED_BLAS -lmwblas

TIMER   = $(SRC_DIR)/timer.c

SRC = $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) 

all: mex

mex: $(SRC) Makefile lbfgsb_wrapper.c
	$(CC) $(FLAGS) $(INCLUDES) lbfgsb_wrapper.c $(SRC)
