CC=clang
CXX=clang++

CFLAGS=-O3
CXXFLAGS=-DSQUIDS_USE_VECTOR_EXTENSIONS=0

CXXFLAGS+=-I${SNOTBUILDPATH}/include/
LDFLAGS+=-L${SNOTBUILDPATH}/lib/

COOLCVMFSROOT=/cvmfs/icecube.opensciencegrid.org/py2-v3_early_access/

CXXFLAGS+=-std=c++11
CXXFLAGS+=-I. -O3
CXXFLAGS+=-I${COOLCVMFSROOT}/RHEL_6_x86_64/include
LDFLAGS+=-L${COOLCVMFSROOT}/RHEL_6_x86_64/lib
LDFLAGS+=-L${COOLCVMFSROOT}/RHEL_6_x86_64/lib64
LDFLAGS+=-lgsl -lgslcblas -lm
LDFLAGS+=-lboost_system -lboost_iostreams -lboost_filesystem -lboost_regex
LDFLAGS+=-lhdf5 -lhdf5_hl
# the following should be found in the include and lib directories
# of the SNOTBUILDPATH
LDFLAGS+=-lNewNuFlux
LDFLAGS+=-lPhysTools
LDFLAGS+=-lSQuIDS -lnuSQuIDS
LDFLAGS+=-lLeptonWeighter
#LDFLAGS+=-lsupc++

CXXFLAGS+=-I../Likelihood
LIKELIHOOD_HEADERS=../Likelihood/autodiff.h ../Likelihood/likelihood.h ../Likelihood/nugenWeighting.h ../Likelihood/particleType.h ../Likelihood/tableio.h ../Likelihood/weighting.h ../Likelihood/lbfgsb/lbfgsb.h ../Likelihood/lbfgsb/interface.h

ICSterile_HELPER_OBJECTS = lbfgsb.o linpack.o analysisWeighting.o dataIO.o Event.o plottingExtras.o
ICSterile_OBJECT = SterileSearchFit.o
ICSterileRemastered_OBJECT = SterileSearchFitRemastered.o

.PHONY: all clean

all: SterileSearchFit SterileSearchFitRemastered

clean:
	rm -rf SterileSearchFit $(ICSterileRemastered_OBJECT) $(ICSterile_HELPER_OBJECTS) $(ICSterile_OBJECT)

lbfgsb.o : ../Likelihood/lbfgsb/lbfgsb.h ../Likelihood/lbfgsb/lbfgsb.c
	$(CC) $(CFLAGS) ../Likelihood/lbfgsb/lbfgsb.c -c -o lbfgsb.o

linpack.o : ../Likelihood/lbfgsb/linpack.c
	$(CC) $(CFLAGS) ../Likelihood/lbfgsb/linpack.c -c -o linpack.o

analysisWeighting.o : analysisWeighting.cpp analysisWeighting.h splinetable.h
	$(CXX) $(CXXFLAGS) analysisWeighting.cpp -c -o analysisWeighting.o

dataIO.o : dataIO.h dataIO.cpp Event.h analysisWeighting.h plottingExtras.h $(LIKELIHOOD_HEADERS)
	$(CXX) $(CXXFLAGS) dataIO.cpp -c -o dataIO.o

Event.o : Event.cpp Event.h
	$(CXX) $(CXXFLAGS) Event.cpp -c -o Event.o

plottingExtras.o : plottingExtras.cpp plottingExtras.h
	$(CXX) $(CXXFLAGS) plottingExtras.cpp -c -o plottingExtras.o

SterileSearchFit.o : SterileSearchFit.cpp $(LIKELIHOOD_HEADERS) analysisWeighting.h Event.h dataIO.h  plottingExtras.h cl_options.h splinetable.h
	$(CXX) $(CXXFLAGS) SterileSearchFit.cpp -c -o SterileSearchFit.o

SterileSearchFitRemastered.o : SterileSearchFitRemastered.cpp $(LIKELIHOOD_HEADERS) analysisWeighting.h Event.h dataIO.h  plottingExtras.h cl_options.h splinetable.h
	$(CXX) $(CXXFLAGS) SterileSearchFitRemastered.cpp -c -o SterileSearchFitRemastered.o

SterileSearchFit : $(ICSterile_HELPER_OBJECTS) $(ICSterile_OBJECT)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lphotospline $(ICSterile_HELPER_OBJECTS) $(ICSterile_OBJECT) -o SterileSearchFit

SterileSearchFitRemastered : $(ICSterile_HELPER_OBJECTS) $(ICSterileRemastered_OBJECT)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lphotospline $(ICSterile_HELPER_OBJECTS) $(ICSterileRemastered_OBJECT) -o SterileSearchFitRemastered

