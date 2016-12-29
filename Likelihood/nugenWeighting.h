#ifndef NUGENWEIGHTING_H
#define NUGENWEIGHTING_H

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "weighting.h"

#include <NewNuFlux/particleType.h>

namespace{ //ugly ODR hack
bool sameGeneration(I3Particle::ParticleType p1, I3Particle::ParticleType p2){
	switch(p1){
		case I3Particle::NuE:
		case I3Particle::NuEBar:
			return(p2==I3Particle::NuE || p2==I3Particle::NuEBar);
		case I3Particle::NuMu:
		case I3Particle::NuMuBar:
			return(p2==I3Particle::NuMu || p2==I3Particle::NuMuBar);
		case I3Particle::NuTau:
		case I3Particle::NuTauBar:
			return(p2==I3Particle::NuTau || p2==I3Particle::NuTauBar);
		default:
			return(false);
	}
}
}

struct nugenSimDetails{
	unsigned long files;
	unsigned long eventsPerFile;
	I3Particle::ParticleType flavor;
	int year;
	double injectionRadius;
	double azimuthMin;
	double azimuthMax;
	double zenithMin;
	double zenithMax;
	double energyMin;
	double energyMax;
	double powerlawIndex;
	double unshadowedFraction; //not actually used for NuGen weighting, but is a property of each sim. set
	
	double generationSolidAngle() const{
		return((azimuthMax-azimuthMin)*(cos(zenithMin)-cos(zenithMax)));
	}
	
	//returns result in cm^2 !
	double generationArea() const{
		return(1e4*boost::math::constants::pi<double>()*injectionRadius*injectionRadius);
	}
	
	bool inBounds(double energy, double zenith, double azimuth/*, double impactParameter*/, I3Particle::ParticleType type, int year){
		return(energy>=energyMin && energy<=energyMax
			   && zenith>=zenithMin && zenith<=zenithMax
			   && azimuth>=azimuthMin && azimuth<azimuthMax
			   //&& impactParameter<=injectionRadius
			   && this->year==year
			   && sameGeneration(flavor,type));
	}
	
	double generationProbability(double energy, double zenith, double azimuth/*, double impactParameter*/, I3Particle::ParticleType type, int year){
		if(!inBounds(energy,zenith,azimuth,/*impactParameter,*/type,year))
			return(0);
		double norm=0;
		if(powerlawIndex!=1)
			norm=(1-powerlawIndex)/(pow(energyMax,1-powerlawIndex)-pow(energyMin,1-powerlawIndex));
		else if(powerlawIndex==1)
			norm=1/log(energyMax/energyMin);
		return(eventsPerFile*files*norm*pow(energy,-powerlawIndex)/(generationArea()*generationSolidAngle()));
	}
};

struct nugenWeighter : public GenericWeighter<nugenWeighter>{
private:
	std::vector<nugenSimDetails> generationSpectra;
public:
	nugenWeighter(){}
	nugenWeighter(std::vector<nugenSimDetails> spectra):generationSpectra(spectra){}
	nugenWeighter(std::initializer_list<nugenSimDetails> spectra):generationSpectra(spectra){}
	
	void addGenerationSpectrum(nugenSimDetails spectrum){ generationSpectra.push_back(spectrum); }
	
	using result_type=double;
	template<typename Event>
	result_type operator()(const Event& e) const{
		assert(!generationSpectra.empty() && "Must specify some generation spectra before computing weights!");
		double denom=0;
		for(auto spectrum : generationSpectra)
			denom+=spectrum.generationProbability(e.primaryEnergy,e.primaryZenith,e.primaryAzimuth,/*e.impactParameter,*/e.primaryType,e.year);
		return(e.interactionProbability/denom);
	}
};

#endif
