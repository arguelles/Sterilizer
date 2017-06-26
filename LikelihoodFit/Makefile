CC=clang
CXX=clang++

CFLAGS=-O3 -fPIC
#CXXFLAGS=-DSQUIDS_USE_VECTOR_EXTENSIONS=0

CXXFLAGS+=-std=c++11 -fPIC
CXXFLAGS+=-I. -O3

COOLCVMFSROOT=/cvmfs/icecube.opensciencegrid.org/py2-v3_early_access/
CXXFLAGS+=-I${COOLCVMFSROOT}/RHEL_6_x86_64/include
LDFLAGS+=-L${COOLCVMFSROOT}/RHEL_6_x86_64/lib
LDFLAGS+=-L${COOLCVMFSROOT}/RHEL_6_x86_64/lib64
CXXFLAGS+=-DSQUIDS_USE_VECTOR_EXTENSIONS=0 -std=c++11

CXXFLAGS+=-I${SNOTBUILDPATH}/include/
LDFLAGS+=-L${SNOTBUILDPATH}/lib/
LDFLAGS+=-lgsl -lgslcblas -lm
LDFLAGS+=-lboost_system -lboost_iostreams -lboost_filesystem -lboost_regex
LDFLAGS+=-lhdf5 -lhdf5_hl
# the following should be found in the include and lib directories
# of the SNOTBUILDPATH
LDFLAGS+=-lNewNuFlux
LDFLAGS+=-lPhysTools
LDFLAGS+=-lSQuIDS -lnuSQuIDS
LDFLAGS+=-lLeptonWeighter
LDFLAGS+=-lsupc++

CXXFLAGS+=-I../Likelihood
LIKELIHOOD_HEADERS=../Likelihood/autodiff.h ../Likelihood/likelihood.h ../Likelihood/nugenWeighting.h ../Likelihood/particleType.h ../Likelihood/tableio.h ../Likelihood/weighting.h ../Likelihood/lbfgsb/lbfgsb.h ../Likelihood/lbfgsb/interface.h

ICSterile_HELPER_OBJECTS = lbfgsb.o linpack.o analysisWeighting.o compactIO.o Event.o SterileSearch.o oversizeWeight.o
ICSterile_OBJECT = SterileSearchTestMinimal.o

.PHONY: all clean

all:  SterileSearchTestMinimal 

clean:
	rm -rf SterileSearchFit $(ICSterileRemastered_OBJECT) $(ICSterile_HELPER_OBJECTS) $(ICSterile_OBJECT)

lbfgsb.o : ../Likelihood/lbfgsb/lbfgsb.h ../Likelihood/lbfgsb/lbfgsb.c
	$(CC) $(CFLAGS) ../Likelihood/lbfgsb/lbfgsb.c -c -o lbfgsb.o

linpack.o : ../Likelihood/lbfgsb/linpack.c
	$(CC) $(CFLAGS) ../Likelihood/lbfgsb/linpack.c -c -o linpack.o

analysisWeighting.o : analysisWeighting.cpp analysisWeighting.h splinetable.h
	$(CXX) $(CXXFLAGS) analysisWeighting.cpp -c -o analysisWeighting.o

compactIO.o : compactIO.h compactIO.cpp Event.h analysisWeighting.h $(LIKELIHOOD_HEADERS)
	$(CXX) $(CXXFLAGS) compactIO.cpp -c -o compactIO.o

Event.o : Event.cpp Event.h
	$(CXX) $(CXXFLAGS) Event.cpp -c -o Event.o

oversizeWeight.o : oversizeWeight.cpp oversizeWeight.h
	$(CXX) $(CXXFLAGS) oversizeWeight.cpp -c -o oversizeWeight.o

SterileSearch.o : SterileSearch.cpp $(LIKELIHOOD_HEADERS) analysisWeighting.h Event.h compactIO.h splinetable.h SterileSearch.h oversizeWeight.h
	$(CXX) $(CXXFLAGS) SterileSearch.cpp -c -o SterileSearch.o

SterileSearchTest.o : SterileSearchTest.cpp $(LIKELIHOOD_HEADERS) analysisWeighting.h Event.h compactIO.h splinetable.h SterileSearch.h
	$(CXX) $(CXXFLAGS) SterileSearchTest.cpp -c -o SterileSearchTest.o

SterileSearchTestMinimal.o : SterileSearchTestMinimal.cpp $(LIKELIHOOD_HEADERS) analysisWeighting.h Event.h compactIO.h splinetable.h SterileSearch.h
	$(CXX) $(CXXFLAGS) SterileSearchTestMinimal.cpp -c -o SterileSearchTestMinimal.o

SterileSearchTestMinimal : $(ICSterile_HELPER_OBJECTS) $(ICSterile_OBJECT)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lphotospline $(ICSterile_HELPER_OBJECTS) $(ICSterile_OBJECT) -o SterileSearchTestMinimal


