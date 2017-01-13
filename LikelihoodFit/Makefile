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

ICSterile_HELPER_OBJECTS = lbfgsb.o linpack.o analysisWeighting.o dataIO.o Event.o SterileSearch.o
ICSterile_OBJECT = SterileSearchDoFit.o

.PHONY: all clean

all: SterileSearchDoFit

clean:
	rm -rf SterileSearchFit $(ICSterileRemastered_OBJECT) $(ICSterile_HELPER_OBJECTS) $(ICSterile_OBJECT)

lbfgsb.o : ../Likelihood/lbfgsb/lbfgsb.h ../Likelihood/lbfgsb/lbfgsb.c
	$(CC) $(CFLAGS) ../Likelihood/lbfgsb/lbfgsb.c -c -o lbfgsb.o

linpack.o : ../Likelihood/lbfgsb/linpack.c
	$(CC) $(CFLAGS) ../Likelihood/lbfgsb/linpack.c -c -o linpack.o

analysisWeighting.o : analysisWeighting.cpp analysisWeighting.h splinetable.h
	$(CXX) $(CXXFLAGS) analysisWeighting.cpp -c -o analysisWeighting.o

dataIO.o : dataIO.h dataIO.cpp Event.h analysisWeighting.h $(LIKELIHOOD_HEADERS)
	$(CXX) $(CXXFLAGS) dataIO.cpp -c -o dataIO.o

Event.o : Event.cpp Event.h
	$(CXX) $(CXXFLAGS) Event.cpp -c -o Event.o

plottingExtras.o : plottingExtras.cpp plottingExtras.h
	$(CXX) $(CXXFLAGS) plottingExtras.cpp -c -o plottingExtras.o

SterileSearch.o : SterileSearch.cpp  $(LIKELIHOOD_HEADERS) analysisWeighting.h Event.h dataIO.h splinetable.h SterileSearch.h
	$(CXX) $(CXXFLAGS) SterileSearch.cpp -c -o SterileSearch.o

SterileSearchDoFit.o : SterileSearchDoFit.cpp $(LIKELIHOOD_HEADERS) analysisWeighting.h Event.h dataIO.h splinetable.h SterileSearch.h
	$(CXX) $(CXXFLAGS) SterileSearchFit.cpp -c -o SterileSearchFit.o

SterileSearchDoFit : $(ICSterile_HELPER_OBJECTS) $(ICSterile_OBJECT)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lphotospline $(ICSterile_HELPER_OBJECTS) $(ICSterile_OBJECT) -o SterileSearchDoFit

