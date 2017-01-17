#include "Event.h"
#include "analysisWeighting.h"

void piecewisePowerlawFlux::read(std::istream& is){
	unsigned int segments;
	double eMin, eMax, norm, power;
	is >> segments;
	if(is.fail())
		throw std::runtime_error("Couldn't read valid number of powerlaw segments");
	for(unsigned int i=0; i<segments; i++){
		is >> eMin >> eMax;
		if(is.fail())
			throw std::runtime_error("Couldn't read segment energy domain");
		unsigned int count;
		is >> count;
		if(is.fail())
			throw std::runtime_error("Couldn't read segment's number of parameters");
		if(count!=2)
			throw std::runtime_error("Segment must have exactly 2 parameters (normalization and index)");
		is >> norm >> power;
		if(is.fail())
			throw std::runtime_error("Couldn't read segment parameters");
		pieces.push_back(powerlaw{eMin,eMax,norm,power});
	}
}

bool operator<(const piecewisePowerlawFlux::powerlaw& p1, const piecewisePowerlawFlux::powerlaw& p2){
	return(p1.eMin<p2.eMin);
}

