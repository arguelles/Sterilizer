#ifndef ANALYSISWEIGHTING_H
#define ANALYSISWEIGHTING_H

#include <cassert>
#include <type_traits>

//#include <NewNuFlux/NewNuFlux.h>
#include <LeptonWeighter/particleType.h>

#include "../Likelihood/autodiff.h"

#include "weighting.h"
#include "splinetable.h"

#define SINGLE_DOMEFF_TEMPLATE

namespace {
	template<class SharedPointer> struct Holder {
		SharedPointer p;
		
		Holder(const SharedPointer &p) : p(p) {}
		Holder(const Holder &other) : p(other.p) {}
		//Holder(Holder &&other) : p(std::move<SharedPointer>(other.p)) {}
		
		void operator () (...) const {}
    };
	
	template<class T> std::shared_ptr<T> to_std_ptr(const boost::shared_ptr<T> &p) {
		typedef Holder<std::shared_ptr<T>> H;
		if(H* h = boost::get_deleter<H, T>(p))
			return h->p;
		return std::shared_ptr<T>(p.get(), Holder<boost::shared_ptr<T>>(p));
	}
	
	template<class T> boost::shared_ptr<T> to_boost_ptr(const std::shared_ptr<T> &p){
		typedef Holder<boost::shared_ptr<T>> H;
		if(H* h = std::get_deleter<H, T>(p))
			return h->p;
		return boost::shared_ptr<T>(p.get(), Holder<std::shared_ptr<T>>(p));
	}
}

template<typename Event, typename DataType>
struct neuaneu : public GenericWeighter<neuaneu<Event,DataType>> {
private:
	DataType ratio;
public:
	neuaneu(DataType ratio): ratio(ratio) {};

	using result_type=DataType;
	result_type operator()(const Event& e) const{
  //std::cout << (int) e.primaryType << std::endl;
  //std::cout << (int) particleType::MuMinus << std::endl;
  if((particleType)e.primaryType == particleType::MuMinus or
     (particleType)e.primaryType == particleType::EMinus  or
     (particleType)e.primaryType == particleType::TauMinus )
		return ratio;
  else
		return 2.-ratio;
	}
};

template<typename Event, typename DataType>
struct AtmosphericZenithVariationCorrectionFactor: public GenericWeighter<AtmosphericZenithVariationCorrectionFactor<Event,DataType>> {
private:
	DataType delta;
	double c0;
  double c1;
  double e0;
  double e1;
public:
	AtmosphericZenithVariationCorrectionFactor(DataType delta): delta(delta),c0(0.4),c1(100.),e0(369.),e1(11279.) {};

	using result_type=DataType;
	result_type operator()(const Event& e) const{
    double ct = cos(e.zenith);
    return 1.0 + (ct+c0)*delta*(1.+(e.energy-e0)/e1)/(1.+exp(-2.*c1*(ct+c0)));
	}
};

struct piecewisePowerlawFlux{
private:
	struct powerlaw{
		double eMin, eMax;
		double norm, power;
		
		bool containsEnergy(double e) const{ return(e>=eMin && e<eMax); }
		double operator()(double e) const{ return(norm*pow(e,power)); }
	};
	friend bool operator<(const powerlaw&, const powerlaw&);
	
	static bool overlap(const powerlaw& p1, const powerlaw& p2){
		if(p1<p2)
			return(p1.eMax<=p2.eMin);
		return(p2.eMax<=p1.eMin);
	}
	
	static bool below_max(double e, const powerlaw& p){ return(e<p.eMax); }
	
	void read(std::istream& is);
	
	std::vector<powerlaw> pieces;
	
public:
	explicit piecewisePowerlawFlux(std::string path){
		std::ifstream file(path.c_str());
		if(!file)
			throw std::runtime_error("Couldn't open "+path+" for reading");
		read(file);
	}
	
	double getFlux(double e) const{
		auto it=std::upper_bound(pieces.begin(),pieces.end(),e,below_max);
		if(it==pieces.end())
			return(0); //treat flux as zero where it is undefined
		if(!it->containsEnergy(e))
			return(0); //treat flux as zero where it is undefined
		return((*it)(e));
	}
};


template<typename T>
struct powerlawWeighter : public GenericWeighter<powerlawWeighter<T>>{
private:
	T index;
public:
	using result_type=T;
	powerlawWeighter(T i):index(i){}
	
	template<typename Event>
	result_type operator()(const Event& e) const{
		return(pow((double)e.primaryEnergy,index));
	}
};

template<typename T, typename Event, typename U>
struct cachedValueWeighter : public GenericWeighter<cachedValueWeighter<T,Event,U>>{
private:
	U Event::* cachedPtr;
public:
	using result_type=T;
	cachedValueWeighter(U Event::* ptr):cachedPtr(ptr){}
	result_type operator()(const Event& e) const{
		return(result_type(e.*cachedPtr));
	}
};

template<typename T, typename Event>
struct FunctionWeighter : public GenericWeighter<FunctionWeighter<T,Event>>{
private:
	std::function<T(const Event&)> func;
public:
	using result_type=T;
	FunctionWeighter(std::function<T(const Event&)> f):func(f){}
	result_type operator()(const Event& e) const{
		return(result_type(func(e)));
	}
};

template<typename T, typename E>
FunctionWeighter<T,E> makeFunctionWeighter(std::function<T(const E&)> f){
	return(FunctionWeighter<T,E>(f));
}


template<typename Event, typename T>
struct boxEnergyFilter : public GenericWeighter<boxEnergyFilter<Event,T>>{
private:
	double min, max;
public:
	using result_type=T;
	boxEnergyFilter(double min, double max):min(min),max(max){}
	
	result_type operator()(const Event& e) const{
		return((e.injectedEnergy>min && e.injectedEnergy<max) ? 1 : 0);
	}
};

//Tilt a spectrum by an incremental powerlaw index about a precomputed median energy
template<typename Event, typename T>
struct powerlawTiltWeighter : public GenericWeighter<powerlawTiltWeighter<Event,T>>{
private:
	double medianEnergy;
	T deltaIndex;
public:
	using result_type=T;
	powerlawTiltWeighter(double me, T dg):
	medianEnergy(me),deltaIndex(dg){}
	
	result_type operator()(const Event& e) const{
		result_type weight=pow(e.injectedEnergy/medianEnergy,-deltaIndex);
		return(weight);
	}
};

namespace DOMEff3{
  using DOMMapType=std::map<unsigned int,std::shared_ptr<Splinetable>>;
	template<typename Event, typename DataType>
	struct domEffWeighter : public GenericWeighter<domEffWeighter<Event,DataType>>{
    private:
      DOMMapType domEffMap;
      DataType logEff;
    public:
      domEffWeighter(DOMMapType domEffMap, DataType deltaEff):
      domEffMap(domEffMap),logEff(log10(0.9*(1.0+deltaEff))){}

      using result_type=DataType;
      result_type operator()(const Event& e) const{
        float cache=e.cachedDOMEff;
        double coordinates[3]={log10(e.energy),cos(e.zenith),logEff};

        auto domcorrection = domEffMap.find(e.year);
        if( domcorrection == domEffMap.end() )
          throw std::runtime_error("Dom efficiency correction for year "+std::to_string(e.year) + " not found.");
        double rate=(*(*domcorrection).second)(coordinates);
        return(pow(10.,rate-cache));
      }
    };

    template<typename Event, int Dim>
    struct domEffWeighter<Event,FD<Dim>> : public GenericWeighter<domEffWeighter<Event,FD<Dim>>>{
    private:
      DOMMapType domEffMap;
      FD<Dim> logEff;
      unsigned int didx;
    public:
      domEffWeighter(DOMMapType domEffMap, FD<Dim> deltaEff):
      domEffMap(domEffMap),
      logEff(log10(0.9*(1.0+deltaEff)))
      {

        const unsigned int n=detail::dimensionExtractor<FD,Dim,double>::nVars(deltaEff);
        for(int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
          if(deltaEff.derivative(i)!=0){
            didx=i;
            break;
          }
        }
      }
      using result_type=FD<Dim>;
      result_type operator()(const Event& e) const{
        float cache=e.cachedDOMEff;
        double rate, derivative;
        double coordinates[3]={log10(e.energy),cos(e.zenith),logEff.value()};

        auto domcorrection = domEffMap.find(e.year);
        if( domcorrection == domEffMap.end() )
          throw std::runtime_error("Dom efficiency correction for year "+std::to_string(e.year) + " not found.");
        rate=(*(*domcorrection).second)(coordinates);
        derivative=((*domcorrection).second)->derivative(coordinates,2);
        derivative*=logEff.derivative(didx);

        result_type r(rate);
        r.setDerivative(derivative,didx);
        return(pow(10.,r-cache));
      }
	};

	//used for initializing the per-event dom-eff related caches
	template<typename Event>
	struct simpleEffRate{
		Splinetable* rate2011;
		double logEff;
		simpleEffRate(Splinetable* r2011, double eff):
		rate2011(r2011),logEff(log10(eff)){}

		void setCache(Event& e) const{
			double coordinates[3]={log10(e.energy),cos(e.zenith),logEff};
			Splinetable* rateTable=nullptr;
			switch(e.year){
				case 2011: rateTable=rate2011; break;
				default: assert(false && "Unexpected year");
			}
			e.cachedDOMEff=(*rateTable)(coordinates);
		}
	};
}

using namespace DOMEff3;


struct DiffuseFitWeighterMaker{
private:
	//median value for dataset 6454 weighted with honda2006_gaisserH3a_elbert_numu
	static constexpr double medianConvEnergy=2020;
	//median value for dataset 6454 weighted with sarcevic_std_gaisserH3a_elbert_numu
	static constexpr double medianPromptEnergy=7887;
  // DOM efficienfy splines;
  using DOMMapType=std::map<unsigned int,std::shared_ptr<Splinetable>>;
  DOMMapType domEffConv_;
  DOMMapType domEffPrompt_;
  bool spline_sets;
public:
	DiffuseFitWeighterMaker(DOMMapType domEffConv_,DOMMapType domEffPrompt_):
    domEffConv_(domEffConv_),domEffPrompt_(domEffPrompt_),spline_sets(true)
	{}

  // default constructor // bad bad
	DiffuseFitWeighterMaker():spline_sets(false){}

  void SetSplines(DOMMapType domEffConv){
    domEffConv_=domEffConv;
    spline_sets=true;
  }

	template<typename DataType>
	std::function<DataType(const Event&)> operator()(const std::vector<DataType>& params) const{
    if(not spline_sets)
      throw std::runtime_error("DFWM object splines not set.");
		assert(params.size()==8);
		//unpack things so we have legible names
		DataType convNorm=params[0];
		DataType promptNorm=params[1];
		DataType astroNorm=params[2];
		DataType CRDeltaGamma=params[3];
		DataType deltaDomEff=params[4];
		DataType piKRatio=params[5];
		DataType NeutrinoAntineutrinoRatio=params[6];
		DataType AtmosphericZenithVariationCorrectionFactorParameter=params[7];

		using cachedWeighter=cachedValueWeighter<DataType,Event,float>;
		//cachedWeighter livetime(&Event::cachedLivetime);
		cachedWeighter convPionFlux(&Event::cachedConvPionWeight);
		cachedWeighter convKaonFlux(&Event::cachedConvKaonWeight);
		cachedWeighter promptFlux(&Event::cachedPromptWeight);

		using domEffW_t = domEffWeighter<Event,DataType>;
		domEffW_t convDOMEff(domEffConv_,deltaDomEff);

		using neuaneu_t = neuaneu<Event,DataType>;
		neuaneu_t neuaneu_w(NeutrinoAntineutrinoRatio);
		using AtmosphericZenithVariationCorrectionFactor_t = AtmosphericZenithVariationCorrectionFactor<Event,DataType>;
		AtmosphericZenithVariationCorrectionFactor_t AtmosphericZenithVariationCorrectionFactor(AtmosphericZenithVariationCorrectionFactorParameter);

		auto conventionalComponent = convNorm*(convPionFlux + piKRatio*convKaonFlux)
		                             *powerlawTiltWeighter<Event,DataType>(medianConvEnergy, CRDeltaGamma)
		                             *convDOMEff*neuaneu_w*AtmosphericZenithVariationCorrectionFactor;

		auto promptComponent = promptNorm*promptFlux
		                       *powerlawTiltWeighter<Event,DataType>(medianPromptEnergy, CRDeltaGamma);

		return (conventionalComponent+promptComponent);
	}
};

#endif
