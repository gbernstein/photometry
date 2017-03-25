# These should come from environment:
# CXX
# OPTFLAGS

# TMV_DIR
# FFTW_DIR
# YAML_DIR
# CFITSIO_DIR
# MKL_DIR (optionally)

# It is assumed that the utilities/ directory is at ../utilities 

export CXX
export OPTFLAGS

# ABS_INCLUDES are absolute paths 
ABS_INCLUDES = -I $(TMV_DIR)/include -I $(YAML_DIR)/include

ifdef MKL_DIR
ABS_INCLUDES += -I $(MKL_DIR)/include
endif


SUBDIRS = 

INCLUDES = -I ../utilities -I ../astrometry

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)

SRC = $(shell ls *.cpp)

OBJ = PhotoMap.o PhotoMapCollection.o SubMap.o PhotoPrior.o PhotoMatch.o PhotoTemplate.o

all:  $(OBJ)

# For building test programs:
UTILITIES := ../utilities
SUBOBJ = $(UTILITIES)/StringStuff.o $(UTILITIES)/Poly2d.o 
LIB_DIRS = -L $(CFITSIO_DIR)/lib -L $(TMV_DIR)/lib -L $(FFTW_DIR)/lib \
	-L $(YAML_DIR)/lib
TMV_LINK := $(shell cat $(TMV_DIR)/share/tmv/tmv-link)
CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)
LIBS = -lm $(LIB_DIRS) -lyaml-cpp -lfftw3 -lcfitsio -ltmv_symband $(TMV_LINK)

testTemplate: testTemplate.o $(OBJ) $(SUBOBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@

###############################################################
## Standard stuff:
###############################################################


subs:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE)); done

depend:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) depend); done
	$(CXX) $(CXXFLAGS) -MM $(SRC) > .$@

clean:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) clean); done
	rm -f *.o *~ *.aux *.log *.dvi core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

export

.PHONY: all install dist depend clean 
