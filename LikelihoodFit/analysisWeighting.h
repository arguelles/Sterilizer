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

struct astroNeutrinoFluxWeighter : public GenericWeighter<astroNeutrinoFluxWeighter>{
private:
	piecewisePowerlawFlux flux;
public:
	astroNeutrinoFluxWeighter(std::string modelPath):flux(modelPath){}
	
	using result_type=double;
	
	template<typename Event>
	result_type operator()(const Event& e) const{
		return(flux.getFlux(e.primaryEnergy));
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

/*template<typename T>
struct powerlawWithExpCutoffWeighter : public GenericWeighter<powerlawWithExpCutoffWeighter<T>>{
private:
	T index;
	T cutoff;
public:
	using result_type=T;
	powerlawWeighter(T i, T c):index(i),cutoff(c){}
	
	template<typename Event>
	result_type operator()(const Event& e) const{
		return(pow((double)e.primaryEnergy,index)*exp(-e.primaryEnergy/cutoff));
	}
};*/

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
	//typename Event::crTiltValues Event::* cachedData;
public:
	using result_type=T;
	powerlawTiltWeighter(double me, T dg/*, typename Event::crTiltValues Event::* c*/):
	medianEnergy(me),deltaIndex(dg)/*,cachedData(c)*/{}
	
	result_type operator()(const Event& e) const{
		//const typename Event::crTiltValues& cache=e.*cachedData;
		result_type weight=pow(e.injectedEnergy/medianEnergy,-deltaIndex);
		return(weight);
	}
};
/*template<typename Event>
struct powerlawTiltWeighter<Event,double> : public GenericWeighter<powerlawTiltWeighter<Event,double>>{
private:
	double medianEnergy;
	double deltaIndex;
	typename Event::crTiltValues Event::* cachedData;
public:
	using result_type=double;
	powerlawTiltWeighter(double me, double dg, typename Event::crTiltValues Event::* c):
	medianEnergy(me),deltaIndex(dg),cachedData(c){}
	
	result_type operator()(const Event& e) const{
		const typename Event::crTiltValues& cache=e.*cachedData;
		if(cache.index==deltaIndex)
			return(cache.weight);
		double weight=pow(e.primaryEnergy/medianEnergy,-deltaIndex);
		cache.index=deltaIndex;
		cache.weight=weight;
		return(weight);
	}
};*/

namespace DOMEff1{

struct DOMEffCorrectionFactory{
	
	struct quadratic{
		double a,b,c;
		
		template<typename T>
		T operator()(T x) const{
			return((a*x+b)*x+c);
		}
	};
	
	template<typename T>
	struct rate{
		T eff;
		T A,Ec,a,b,c,d,f;
		
		rate():A(0.0){}
		rate(T eff, T A, T Ec, T a, T b, T c, T d, T f):
		eff(eff),A(A),Ec(Ec),a(a),b(b),c(c),d(d),f(f){}
		
		T operator()(T E) const{
			T Es=E/Ec;
			T lE=log10(E);
			return(A/(pow(Es,((a*lE+b)*lE+c)*lE+d)+pow(Es,f)));
		}
	};
	
	quadratic A,Ec,a,b,c,d,f;
	
	template<typename T>
	rate<T> rateForEfficiency(T e) const{
		std::cout << "rate for efficiency " << e << ":\n";
		std::cout << ' ' << A(e) << '\n';
		std::cout << ' ' << Ec(e) << '\n';
		std::cout << ' ' << a(e) << '\n';
		std::cout << ' ' << b(e) << '\n';
		std::cout << ' ' << c(e) << '\n';
		std::cout << ' ' << d(e) << '\n';
		std::cout << ' ' << f(e) << '\n';
		return(rate<T>{e,A(e),Ec(e),a(e),b(e),c(e),d(e),f(e)});
	}
};

	extern DOMEffCorrectionFactory convDOMEffParams;
	extern DOMEffCorrectionFactory promptDOMEffParams;
	extern DOMEffCorrectionFactory astroDOMEffParams;
	
	template<typename Event, typename DataType>
	struct domEffWeighter : public GenericWeighter<domEffWeighter<Event,DataType>>{
	private:
		DOMEffCorrectionFactory::rate<DataType> rate;
		typename Event::domEffValues Event::* cachedData;
	public:
		domEffWeighter(DOMEffCorrectionFactory::rate<DataType> r, typename Event::domEffValues Event::* c):
		rate(r),cachedData(c){}
		
		using result_type=DataType;
		result_type operator()(const Event& e) const{
			const typename Event::domEffValues& cache=e.*cachedData;
			result_type newRate=rate(result_type(e.energy));
			return(newRate/cache.baseRate);
		}
	};
	/*template<typename Event>
	struct domEffWeighter<Event,double> : public GenericWeighter<domEffWeighter<Event,double>>{
	private:
		DOMEffCorrectionFactory::rate<double> rate;
		typename Event::domEffValues Event::* cachedData;
	public:
		domEffWeighter(DOMEffCorrectionFactory::rate<double> r, typename Event::domEffValues Event::* c):
		rate(r),cachedData(c){}
		
		using result_type=double;
		result_type operator()(const Event& e) const{
			const typename Event::domEffValues& cache=e.*cachedData;
			if(cache.lastEff==rate.eff)
				return(cache.lastRate/cache.baseRate);
			double newRate=rate(e.energy);
			cache.lastEff=rate.eff;
			cache.lastRate=newRate;
			return(newRate/cache.baseRate);
		}
	};*/
}

namespace DOMEff2{
	
	struct DOMEffCorrectionFactory{
		
		struct quadratic{
			double a,b,c;
			
			template<typename T>
			T operator()(T x) const{
				return((a*x+b)*x+c);
			}
			
			static quadratic fit(double x1, double y1, double x2, double y2, double x3, double y3){
				double a=((y3-y1)*(x2-x1)-(y2-y1)*(x3-x1))/((x3*x3-x1*x1)*(x2-x1)-(x2*x2-x1*x1)*(x3-x1));
				double b=(y2-y1+a*(x1*x1-x2*x2))/(x2-x1);
				double c=y1-a*x1*x1-b*x1;
				return(quadratic{a,b,c});
			}
		};
		
		template<typename T>
		struct rate{
			T eff;
			T A,Ec,a,b,c,d,f;
			T B,bm,bs;
			T C,cm,cs;
			T D,dm,ds;
			bool normalCorrections;
			
			rate():A(0.0){}
			rate(T eff, T A, T Ec, T a, T b, T c, T d, T f,
				 T B, T bm, T bs,
				 T C, T cm, T cs,
				 T D, T dm, T ds):
			eff(eff),A(A),Ec(Ec),a(a),b(b),c(c),d(d),f(f),
			B(B),bm(bm),bs(bs),
			C(C),cm(cm),cs(cs),
			D(D),dm(dm),ds(ds),
			normalCorrections(B!=0 || C!=0 || D!=0){}
			
			static T normal(T m, T s, T x){
				T z=(log10(x)-m)/s;
				return(exp(-z*z)/(s*x));
			}
			
			T operator()(T E) const{
				T Es=E/Ec;
				T lE=log10(E);
				T r=A/(pow(Es,((a*lE+b)*lE+c)*lE+d)+pow(Es,f));
				if(normalCorrections)
					r+=B*normal(bm,bs,E)
					  +C*normal(cm,cs,E)
					  +D*normal(dm,ds,E);
				return(r);
			}
		};
		
		quadratic A,Ec,a,b,c,d,f;
		quadratic B,bm,bs,C,cm,cs,D,dm,ds;
		
		template<typename T>
		rate<T> rateForEfficiency(T e) const{
			/*std::cout << "rate for efficiency " << e << ":\n";
			std::cout << ' ' << A(e) << '\n';
			std::cout << ' ' << Ec(e) << '\n';
			std::cout << ' ' << a(e) << '\n';
			std::cout << ' ' << b(e) << '\n';
			std::cout << ' ' << c(e) << '\n';
			std::cout << ' ' << d(e) << '\n';
			std::cout << ' ' << f(e) << '\n';
			std::cout << ' ' << B(e) << '\n';
			std::cout << ' ' << bm(e) << '\n';
			std::cout << ' ' << bs(e) << '\n';
			std::cout << ' ' << C(e) << '\n';
			std::cout << ' ' << cm(e) << '\n';
			std::cout << ' ' << cs(e) << '\n';
			std::cout << ' ' << D(e) << '\n';
			std::cout << ' ' << dm(e) << '\n';
			std::cout << ' ' << ds(e) << '\n';*/
			return(rate<T>{e,A(e),Ec(e),a(e),b(e),c(e),d(e),f(e),B(e),bm(e),bs(e),C(e),cm(e),cs(e),D(e),dm(e),ds(e)});
		}
	};
	
	//icky, but hard-code dates for speed
	template<typename DataType>
	struct DOMEffRate{
	private:
		DOMEffCorrectionFactory::rate<DataType> rate2010;
		DOMEffCorrectionFactory::rate<DataType> rate2011;
	public:
		const DataType efficiency;
		
		DOMEffRate(DOMEffCorrectionFactory::rate<DataType> r2010,
				   DOMEffCorrectionFactory::rate<DataType> r2011):
		rate2010(r2010),rate2011(r2011),efficiency(rate2010.eff){
			assert(rate2010.eff==rate2011.eff);
		}
		
		DOMEffRate(DOMEffCorrectionFactory::rate<DataType> r2011):
		rate2011(r2011),efficiency(rate2011.eff){}

		DataType rate(DataType energy, unsigned int year) const{
			switch(year){
				case 2010:
					return(rate2010(energy));
				case 2011:
					return(rate2010(energy));
				default:
					throw std::runtime_error("Unsupported year for DOM efficiency correction");
			}
		}
	};
	
	//extern DOMEffCorrectionFactory convDOMEffParams2010;
	//extern DOMEffCorrectionFactory promptDOMEffParams2010;
	//extern DOMEffCorrectionFactory astroDOMEffParams2010;
	extern DOMEffCorrectionFactory convDOMEffParams2011;
	extern DOMEffCorrectionFactory promptDOMEffParams2011;
	extern DOMEffCorrectionFactory astroDOMEffParams2011;
	
	template<typename Event, typename DataType>
	struct domEffWeighter : public GenericWeighter<domEffWeighter<Event,DataType>>{
	private:
		DOMEffRate<DataType> rate;
		typename Event::domEffValues Event::* cachedData;
	public:
		domEffWeighter(DOMEffRate<DataType> r, typename Event::domEffValues Event::* c):
		rate(r),cachedData(c){}
		
		using result_type=DataType;
		result_type operator()(const Event& e) const{
			const typename Event::domEffValues& cache=e.*cachedData;
			return(rate.rate(e.energy,e.year)/cache.baseRate);
		}
	};
	/*template<typename Event>
	struct domEffWeighter<Event,double> : public GenericWeighter<domEffWeighter<Event,double>>{
	private:
		DOMEffRate<double> rate;
		typename Event::domEffValues Event::* cachedData;
	public:
		domEffWeighter(DOMEffRate<double> r, typename Event::domEffValues Event::* c):
		rate(r),cachedData(c){}
		
		using result_type=double;
		result_type operator()(const Event& e) const{
			const typename Event::domEffValues& cache=e.*cachedData;
			if(cache.lastEff==rate.efficiency)
				return(cache.lastRate/cache.baseRate);
			double newRate=rate.rate(e.energy,e.year);
			cache.lastEff=rate.efficiency;
			//TODO: why not store rate with the base divided out in cache?
			cache.lastRate=newRate;
			return(newRate/cache.baseRate);
		}
	};*/
}

namespace DOMEff3{
	template<typename Event, typename DataType>
	struct domEffWeighter : public GenericWeighter<domEffWeighter<Event,DataType>>{
	private:
		//Splinetable* rate2010;
		Splinetable* rate2011;
		DataType logEff;
		typename Event::domEffValues Event::* cachedData;
	public:
		//domEffWeighter(Splinetable* r2010, Splinetable* r2011, DataType deltaEff, typename Event::domEffValues Event::* c):
		//rate2010(r2010),rate2011(r2011),logEff(log10(0.9*(1.0+deltaEff))),cachedData(c){}

		domEffWeighter(Splinetable* r2011, DataType deltaEff, typename Event::domEffValues Event::* c):
		rate2011(r2011),logEff(log10(0.9*(1.0+deltaEff))),cachedData(c){}

		using result_type=DataType;
		result_type operator()(const Event& e) const{
			const typename Event::domEffValues& cache=e.*cachedData;
			double rate;
			double coordinates[3]={log10(e.energy),cos(e.zenith),logEff};
			switch(e.year){
		//		case 2010: rate=(*rate2010)(coordinates); break;
				case 2011: rate=(*rate2011)(coordinates); break;
			}
			//return(rate/cache.baseRate);
			return(pow(10.,rate-cache.baseRate));
		}
	};

	template<typename Event, int Dim>
	struct domEffWeighter<Event,FD<Dim>> : public GenericWeighter<domEffWeighter<Event,FD<Dim>>>{
	private:
		Splinetable* rate2011;
		FD<Dim> logEff;
		unsigned int didx;
		typename Event::domEffValues Event::* cachedData;
	public:
    /*
		domEffWeighter(Splinetable* r2010, Splinetable* r2011, FD<Dim> deltaEff, typename Event::domEffValues Event::* c):
		rate2010(r2010),rate2011(r2011),
		logEff(log10(0.9*(1.0+deltaEff))),
		cachedData(c){

			const unsigned int n=detail::dimensionExtractor<FD,Dim,double>::nVars(deltaEff);
			for(int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
				if(deltaEff.derivative(i)!=0){
					didx=i;
					break;
				}
			}
		}*/

		domEffWeighter(Splinetable* r2011, FD<Dim> deltaEff, typename Event::domEffValues Event::* c):
		rate2011(r2011),
		logEff(log10(0.9*(1.0+deltaEff))),
		cachedData(c){

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
			const typename Event::domEffValues& cache=e.*cachedData;
			double rate, derivative;
			double coordinates[3]={log10(e.energy),cos(e.zenith),logEff.value()};
			switch(e.year){
		//		case 2010:
		//			rate=(*rate2010)(coordinates);
		//			derivative=rate2010->derivative(coordinates,2);
		//			//derivative=rate2010->derivativeFast(coordinates,2);
		//			break;
				case 2011:
					rate=(*rate2011)(coordinates);
					derivative=rate2011->derivative(coordinates,2);
					//derivative=rate2011->derivativeFast(coordinates,2);
					break;
			}
			//std::cout << ' ' << rate << ' ' << derivative << " (" << cache.baseRate << ')' << std::endl;
			derivative*=logEff.derivative(didx);
			result_type r(rate);
			r.setDerivative(derivative,didx);
			//std::cout << ' ' << r << std::endl;
			//return(rate/cache.baseRate);
			return(pow(10.,r-cache.baseRate));
		}
	};

	//extern std::unique_ptr<Splinetable> domEffConv2010;

	extern std::unique_ptr<Splinetable> domEffConv2011;
	extern std::unique_ptr<Splinetable> domEffConvPion2011;
	extern std::unique_ptr<Splinetable> domEffConvKaon2011;
	extern std::unique_ptr<Splinetable> domEffPrompt2011;

	//std::unique_ptr<Splinetable> domEffPrompt2010(new Splinetable("/home/carguelles/work/IC_sterile/LikelihoodFit/dom_eff_fits/prompt_IC79.fits"));
	//std::unique_ptr<Splinetable> domEffPrompt2011(new Splinetable("/home/carguelles/work/IC_sterile/LikelihoodFit/dom_eff_fits/prompt_IC86.fits"));
	//std::unique_ptr<Splinetable> domEffAstro2010(new Splinetable("/home/carguelles/work/IC_sterile/LikelihoodFit/dom_eff_fits/astro_IC79.fits"));
	//std::unique_ptr<Splinetable> domEffAstro2011(new Splinetable("/home/carguelles/work/IC_sterile/LikelihoodFit/dom_eff_fits/astro_IC86.fits"));

	//used for initializing the per-event dom-eff related caches
	template<typename Event>
	struct simpleEffRate{
		//Splinetable* rate2010;
		Splinetable* rate2011;
		double logEff;
		typename Event::domEffValues Event::* cachedData;

		simpleEffRate(Splinetable* r2011, double eff, typename Event::domEffValues Event::* c):
		rate2011(r2011),logEff(log10(eff)),cachedData(c){}

		//simpleEffRate(Splinetable* r2010, Splinetable* r2011, double eff, typename Event::domEffValues Event::* c):
		//rate2010(r2010),rate2011(r2011),logEff(log10(eff)),cachedData(c){}
		//rate2010(r2010),rate2011(r2011),logEff(log10(0.9*(1.0+deltaEff))),cachedData(c){}

		void setCache(Event& e) const{
			typename Event::domEffValues& cache=e.*cachedData;
			double coordinates[3]={log10(e.energy),cos(e.zenith),logEff};
			Splinetable* rateTable=nullptr;
			switch(e.year){
		//		case 2010: rateTable=rate2010; break;
				case 2011: rateTable=rate2011; break;
				default: assert(false && "Unexpected year");
			}
			cache.baseRate=(*rateTable)(coordinates);
		}
	};
}

/*
using namespace DOMEff3;

namespace DOMEff3{
	template<typename Event, typename DataType>
	struct domEffWeighter : public GenericWeighter<domEffWeighter<Event,DataType>>{
	private:
		Splinetable* rate2010;
		Splinetable* rate2011;
		DataType logEff;
		typename Event::domEffValues Event::* cachedData;
	public:
		domEffWeighter(Splinetable* r2010, Splinetable* r2011, DataType deltaEff, typename Event::domEffValues Event::* c):
		rate2010(r2010),rate2011(r2011),logEff(log10(0.9*(1+deltaEff))),cachedData(c){}
		
		using result_type=DataType;
		result_type operator()(const Event& e) const{
			const typename Event::domEffValues& cache=e.*cachedData;
			double rate;
			double coordinates[2]={log10(e.energy),logEff};
			switch(e.year){
				case 2010: rate=(*rate2010)(coordinates); break;
				case 2011: rate=(*rate2011)(coordinates); break;
			}
			//return(rate/cache.baseRate);
			return(pow(10.,rate-cache.baseRate));
		}
	};
	
	template<typename Event, int Dim>
	struct domEffWeighter<Event,FD<Dim>> : public GenericWeighter<domEffWeighter<Event,FD<Dim>>>{
	private:
		Splinetable* rate2010;
		Splinetable* rate2011;
		FD<Dim> logEff;
		unsigned int didx;
		typename Event::domEffValues Event::* cachedData;
	public:
		domEffWeighter(Splinetable* r2010, Splinetable* r2011, FD<Dim> deltaEff, typename Event::domEffValues Event::* c):
		rate2010(r2010),rate2011(r2011),
		logEff(log10(0.9*(1+deltaEff))),
		cachedData(c){
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
			const typename Event::domEffValues& cache=e.*cachedData;
			double rate, derivative;
			double coordinates[2]={log10(e.energy),logEff.value()};
			switch(e.year){
				case 2010:
					rate=(*rate2010)(coordinates);
					derivative=rate2010->derivative(coordinates,1);
					break;
				case 2011:
					rate=(*rate2011)(coordinates);
					derivative=rate2011->derivative(coordinates,1);
					break;
			}
			//std::cout << ' ' << rate << ' ' << derivative << " (" << cache.baseRate << ')' << std::endl;
			derivative*=logEff.derivative(didx);
			result_type r(rate);
			r.setDerivative(derivative,didx);
			//std::cout << ' ' << r << std::endl;
			//return(rate/cache.baseRate);
			return(pow(10.,r-cache.baseRate));
		}
	};
	
	extern std::unique_ptr<Splinetable> domEffConv2010;
	extern std::unique_ptr<Splinetable> domEffConv2011;
	extern std::unique_ptr<Splinetable> domEffPrompt2010;
	extern std::unique_ptr<Splinetable> domEffPrompt2011;
	extern std::unique_ptr<Splinetable> domEffAstro2010;
	extern std::unique_ptr<Splinetable> domEffAstro2011;
	
	//used for initializing the per-event dom-eff related caches
	template<typename Event>
	struct simpleEffRate{
		Splinetable* rate2010;
		Splinetable* rate2011;
		double logEff;
		typename Event::domEffValues Event::* cachedData;
		
		simpleEffRate(Splinetable* r2010, Splinetable* r2011, double eff, typename Event::domEffValues Event::* c):
		rate2010(r2010),rate2011(r2011),logEff(log10(eff)),cachedData(c){}
		
		void setCache(Event& e) const{
			typename Event::domEffValues& cache=e.*cachedData;
			double coordinates[2]={log10(e.energy),logEff};
			Splinetable* rateTable=nullptr;
			switch(e.year){
				case 2010: rateTable=rate2010; break;
				case 2011: rateTable=rate2011; break;
				default: assert(false && "Unexpected year");
			}
			cache.baseRate=(*rateTable)(coordinates);
		}
	};
}
*/

using namespace DOMEff3;

template<typename Event, typename T>
struct approximatePowerlawDOMEffCorrector : public GenericWeighter<approximatePowerlawDOMEffCorrector<Event,T>>{
private:
	T deltaIndex, deltaEff;
public:
	using result_type=T;
	approximatePowerlawDOMEffCorrector(T dg, T de):
	deltaIndex(dg),deltaEff(de){}
	
	result_type operator()(const Event& e) const{
		result_type weight=pow(1+deltaEff,deltaIndex);
		return(weight);
	}
};

struct DiffuseFitWeighterMaker{
private:
	//median value for dataset 6454 weighted with honda2006_gaisserH3a_elbert_numu
	static constexpr double medianConvEnergy=2020;
	//median value for dataset 6454 weighted with sarcevic_std_gaisserH3a_elbert_numu
	static constexpr double medianPromptEnergy=7887;
	//value for dataset 6454 weighted with honda2006_gaisserH3a_elbert_numu using MuEx fit from [3e2,1e5]
	/*static constexpr double convThresholdEnergy=240;
	//value for dataset 6454 weighted with honda2006_gaisserH3a_elbert_numu using MuEx fit from [3e2,1e5]
	static constexpr double promptThresholdEnergy=200;
	//value for dataset 6454 weighted with E^-2 using MuEx fit from [3e2,1e6]
	static constexpr double astroThresholdEnergy=195;*/

	//limiting values used for sensitivity calculation
	double minAstroEnergy, maxAstroEnergy;

public:
	DiffuseFitWeighterMaker():
	minAstroEnergy(-std::numeric_limits<double>::infinity()),
	maxAstroEnergy(std::numeric_limits<double>::infinity())
	{}

	DiffuseFitWeighterMaker(double minAstro, double maxAstro):
	minAstroEnergy(minAstro),
	maxAstroEnergy(maxAstro)
	{}

	double getMinAstroEnergy(){ return(minAstroEnergy); }
	void setMinAstroEnergy(double me){ minAstroEnergy=me; }
	double getMaxAstroEnergy(){ return(maxAstroEnergy); }
	void setMaxAstroEnergy(double me){ maxAstroEnergy=me; }

	template<typename DataType>
	std::function<DataType(const Event&)> operator()(const std::vector<DataType>& params) const{
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

		using cachedWeighter=cachedValueWeighter<DataType,Event,double>;
		//cachedWeighter livetime(&Event::cachedLivetime);
		cachedWeighter convPionFlux(&Event::cachedConvPionWeight);
		cachedWeighter convKaonFlux(&Event::cachedConvKaonWeight);
		cachedWeighter promptFlux(&Event::cachedPromptWeight);
		cachedWeighter astroFlux(&Event::cachedAstroWeight);

		using domEffW_t = domEffWeighter<Event,DataType>;
#ifdef SINGLE_DOMEFF_TEMPLATE
		domEffW_t convDOMEff(domEffConv2011.get(),deltaDomEff,&Event::cachedConvDOMEff);
		domEffW_t promptDOMEff(domEffPrompt2011.get(),deltaDomEff,&Event::cachedPromptDOMEff);
#else
		domEffW_t convKaonDOMEff(domEffConvKaon2011.get(),deltaDomEff,&Event::cachedConvDOMEff);
		domEffW_t convPionDOMEff(domEffConvPion2011.get(),deltaDomEff,&Event::cachedConvDOMEff);
		domEffW_t promptDOMEff(domEffPrompt2011.get(),deltaDomEff,&Event::cachedPromptDOMEff);
#endif

		//domEffW_t convDOMEff(domEffConv2010.get(),domEffConv2011.get(),deltaDomEff,&Event::cachedConvDOMEff);
		//domEffW_t astroDOMEff(domEffAstro2010.get(),domEffAstro2011.get(),deltaDomEff,&Event::cachedAstroDOMEff);
		
		using neuaneu_t = neuaneu<Event,DataType>;
		neuaneu_t neuaneu_w(NeutrinoAntineutrinoRatio);
		using AtmosphericZenithVariationCorrectionFactor_t = AtmosphericZenithVariationCorrectionFactor<Event,DataType>;
		AtmosphericZenithVariationCorrectionFactor_t AtmosphericZenithVariationCorrectionFactor(AtmosphericZenithVariationCorrectionFactorParameter);
#ifdef SINGLE_DOMEFF_TEMPLATE
		auto conventionalComponent = convNorm*(convPionFlux + piKRatio*convKaonFlux)
		                             *powerlawTiltWeighter<Event,DataType>(medianConvEnergy, CRDeltaGamma)
		                             *convDOMEff*neuaneu_w*AtmosphericZenithVariationCorrectionFactor;
#else
		auto conventionalComponent = convNorm*(convPionFlux*convPionDOMEff + piKRatio*convKaonFlux*convKaonDOMEff)
		                             *powerlawTiltWeighter<Event,DataType>(medianConvEnergy, CRDeltaGamma)
		                             *neuaneu*AtmosphericZenithVariationCorrectionFactor;
#endif

		//auto promptComponent = promptNorm*promptFlux
//                       *powerlawTiltWeighter<Event,DataType>(medianPromptEnergy, CRDeltaGamma)
		//                       *promptDOMEff;

		auto promptComponent = promptNorm*promptFlux
		                       *powerlawTiltWeighter<Event,DataType>(medianPromptEnergy, CRDeltaGamma);

    /*
		auto astroComponent = astroNorm*astroFlux
		                      *powerlawTiltWeighter<Event,DataType>(1e5, astroDeltaGamma)*approximatePowerlawDOMEffCorrector<Event,DataType>(astroDeltaGamma,deltaDomEff)
		                      *boxEnergyFilter<Event,DataType>(minAstroEnergy,maxAstroEnergy)
		                      *astroDOMEff;
		
		return(conventionalComponent+promptComponent+astroComponent);
    */

		return (conventionalComponent+promptComponent);
	}
};

#endif
