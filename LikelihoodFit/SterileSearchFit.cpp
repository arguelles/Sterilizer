#include <iterator>
#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <chrono>
#include <queue>

#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <LeptonWeighter/lepton_weighter.h>
#include <LeptonWeighter/event.h>
#include <LeptonWeighter/particleType.h>

#include <PhysTools/gnuplot.h>
#include <PhysTools/plottable_histograms.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>

#define NO_STD_OUTPUT
#define COG_CUT

#include "cl_options.h"
#include "likelihood.h"
#include "plottingExtras.h"
#include "Event.h"
#include "analysisWeighting.h"
#include "dataIO.h"
#include "runspec.h"
#include "oversizeWeight.h"

//global settings
size_t evalThreads=1;
bool quiet=false;

bool sameGeneration(particleType p1, particleType p2){
	switch(p1){
		case particleType::NuE:
		case particleType::NuEBar:
			return(p2==particleType::NuE || p2==particleType::NuEBar);
		case particleType::NuMu:
		case particleType::NuMuBar:
			return(p2==particleType::NuMu || p2==particleType::NuMuBar);
		case particleType::NuTau:
		case particleType::NuTauBar:
			return(p2==particleType::NuTau || p2==particleType::NuTauBar);
		default:
			return(false);
	}
}

//Find the median of a weighted set the stupid way: trial and error
template<typename Container, typename Weighter>
double findMedianEnergy(const Container& data, Weighter w, float Event::* energy){
	double min=std::numeric_limits<double>::max(), max=-std::numeric_limits<double>::max();
	double totalWeight=0;
	for(auto& e : data){
		if(e.*energy<min)
			min=e.*energy;
		if(e.*energy>max)
			max=e.*energy;
		totalWeight+=w(e);
	}
	//std::cout << "total weight/rate is " << totalWeight << std::endl;
	while((max-min)>.0001*min){
		double guess=(min+max)/2, weight=0;
		//std::cout << " trying " << guess;
		for(auto& e : data){
			if(e.*energy<guess)
				weight+=w(e);
		}
		//std::cout << ": rate is " << weight << std::endl;
		if(weight>totalWeight/2)
			max=guess;
		else
			min=guess;
	}
	return((min+max)/2);
}

template<typename ContainerType, typename HistType, typename BinnerType>
void bin(const ContainerType& data, HistType& hist, const BinnerType& binner){
	for(const Event& event : data)
		binner(hist,event);
}

//This doesn't belong in here, it is used only for plots
void normalizeSlices(phys_tools::histograms::histogram<2>& h){
	for(unsigned int y=0; y<h.getBinCount(1); y++){
		double total=0.0;
		for(unsigned int x=0; x<h.getBinCount(0); x++)
			total+=h(x,y);
		if(total){
			for(unsigned int x=0; x<h.getBinCount(0); x++)
				h(x,y)/=total;
		}
	}
}

struct fitResult{
	std::vector<double> params;
	double likelihood;
	unsigned int nEval, nGrad;
	bool succeeded;
};

template<typename T>
double FGSL(const gsl_vector *v, void *params){
     T* p = (T*) params;
     return (p->F)(v);
}

template<typename T>
void DFGSL(const gsl_vector *v, void *params, gsl_vector *df){
     T* p = (T*) params;
     (p->DF)(v,df);
}

template<typename T>
void FDFGSL(const gsl_vector *v, void *params, double *f, gsl_vector *df){
     T* p = (T*) params;
     (p->FDF)(v,f,df);
}

template<typename FuncType>
class GSL_Function {
private:
	FuncType func;
public:
	GSL_Function(FuncType f):func(f){}

	virtual double F(const gsl_vector *v){
	    std::vector<double> x;
	    for(unsigned int i = 0; i < v->size; i++){
	    	x.push_back(gsl_vector_get(v,i));
	    }
	    return -(func.template evaluateLikelihood<double>(x));
	}

	virtual void DF(const gsl_vector *v, gsl_vector *df){
	    const size_t size=v->size;

	    std::vector<FD<FuncType::DerivativeDimension>> x(size);
	    for(size_t i=0; i<size; i++)
		    x[i]=FD<FuncType::DerivativeDimension>(gsl_vector_get(v,i),i);
	    FD<FuncType::DerivativeDimension> result=func.template evaluateLikelihood<FD<FuncType::DerivativeDimension>>(x);

	    for(unsigned int i=0; i<size; i++)
		    gsl_vector_set(df,i,-result.derivative(i)); //note negation!
	}

	void FDF(const gsl_vector *v, double *f, gsl_vector *df)
	{
	     *f = F(v);
	     DF(v, df);
	}
};

template<typename LikelihoodType>
fitResult doFitGSL(LikelihoodType& likelihood, const std::vector<double>& seed,
	      std::vector<unsigned int> indicesToFix={}){

    using namespace likelihood;

    int size = seed.size();
    size_t iter = 0;
    size_t iter_max = 500;
    int status;

    #ifdef _USE_DGSL_
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    gsl_multimin_function_fdf func;
    #else
    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function func;
    #endif

    auto gsl_f = GSL_Function<LikelihoodType>(likelihood);

    func.n = size;
    func.f = &FGSL<GSL_Function<LikelihoodType>>;
    #ifdef _USE_DGSL_
    func.df = &DFGSL<GSL_Function<LikelihoodType>>;
    func.fdf = &FDFGSL<GSL_Function<LikelihoodType>>;
    #endif
    func.params = (void *) &gsl_f;

    gsl_vector *x = gsl_vector_alloc (size);
    for(unsigned int i = 0; i < size; i ++)
      gsl_vector_set (x, i, seed[i]);

    #ifdef _USE_DGSL_
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, size);

    gsl_multimin_fdfminimizer_set (s, &func, x, 0.01, 1e-12);
    #else
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, size);
    gsl_vector *xini = gsl_vector_alloc(size);
    for(unsigned int i = 0; i < size; i ++)
      gsl_vector_set (xini, i, 0.005);
    for(unsigned int i : indicesToFix)
      gsl_vector_set (xini, i, 0.0);

    gsl_multimin_fminimizer_set (s, &func, x, xini);
    #endif

    double fval = 0;
    do
    {
        iter++;
        #ifdef _USE_DGSL_
        status = gsl_multimin_fdfminimizer_iterate (s);
        status = gsl_multimin_test_gradient (s->gradient, 1e-12);
        #else
        gsl_multimin_fminimizer_iterate (s);
        //if ( abs(fval - s->fval) < 0.001 ){
        if ( gsl_multimin_fminimizer_size(s)< 1.0e-5 ){
            //std::cout << fval << " " << s->fval << " " << gsl_multimin_fminimizer_size(s) << std::endl;
            status = GSL_SUCCESS;
        }
        else {
            fval = s->fval;
            status = GSL_CONTINUE;
        }
        #endif

        //if (status)
        //break;

        // monitoring stuff DEBUG
        #ifdef _USE_DGSL_
        //std::cout << "status " << status << " " << iter << " ";
        for( int i = 0; i < size; i ++)
            std::cout << gsl_vector_get(s->gradient, i) << " ";
        //std::cout << s->f << std::endl;
        #else
        std::cout << "status " << status << " " << iter << " ";
        std::cout << s->fval << " " << gsl_multimin_fminimizer_size(s) << std::endl;
        #endif

        /*
        if (status != GSL_CONTINUE){
        std::cout << "Stop at:" << std::endl;
        std::cout << iter << " ";
        for( int i = 0; i < size; i ++)
            std::cout << gsl_vector_get(s->x, i) << " ";
        std::cout << s->f << std::endl;
        }
        */
    }
    while (status == GSL_CONTINUE && iter < iter_max);

    fitResult result;
    result.succeeded=(bool)GSL_SUCCESS;
    #ifdef _USE_DGSL_
    result.likelihood=s->f;
    #else
    result.likelihood=s->fval;
    #endif
    std::vector<double> rminloc(size);
    for(int i = 0; i < size; i++)
      rminloc[i] = gsl_vector_get(s->x,i);
    //std::vector<double> minloc = {rminloc[0],0,0,rminloc[1],rminloc[2],rminloc[3]};
    result.params=rminloc;//minloc;
    #ifdef _USE_DGSL_
    gsl_multimin_fdfminimizer_free (s);
    #else
    gsl_multimin_fminimizer_free (s);
    #endif
    gsl_vector_free (x);

    return(result);
}
///Maximize a likelihood using the LBFGSB minimization algorithm
///\param likelihood The likelihood to maximize
///\param seed The seed values for all likelihood parameters
///\param indicesToFix The indices of likelihood parameters to hold constant at their seed values
template<typename LikelihoodType>
fitResult doFitLBFGSB(LikelihoodType& likelihood, const std::vector<double>& seed,
					  std::vector<unsigned int> indicesToFix={}){
	using namespace likelihood;
	
	LBFGSB_Driver minimizer;
	minimizer.addParameter(seed[0],.001,0.0);
	minimizer.addParameter(seed[1],.001,0.0);
	minimizer.addParameter(seed[2],.01,0.0);
	minimizer.addParameter(seed[3],.005);
	minimizer.addParameter(seed[4],.005,-.1,.3/*,.039,.159*/);
	minimizer.addParameter(seed[5],.01,0.0);
	//minimizer.addParameter(seed[6],.1,0.0,2.0);
	minimizer.addParameter(seed[6],.001,0.0,2.0);
	minimizer.addParameter(seed[7],.001,-1.0,1.0);
	
	for(auto idx : indicesToFix)
		minimizer.fixParameter(idx);
	
	//minimizer.setChangeTolerance(1e-4);
	minimizer.setChangeTolerance(1e-5);
	//minimizer.setChangeTolerance(0);
  //minimizer.setGradientTolerance(1.0e-15);
	minimizer.setHistorySize(20);
	
	fitResult result;
	result.succeeded=minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));
	result.likelihood=minimizer.minimumValue();
	result.params=minimizer.minimumPosition();
	result.nEval=minimizer.numberOfEvaluations();
	result.nGrad=minimizer.numberOfEvaluations(); //gradient is always eval'ed with function
		
	return(result);
}

///Wraps a vector of doubles for convenient printing
struct prettyPrint{
private:
	std::vector<double>& v;
public:
	prettyPrint(std::vector<double>& v):v(v){}
	
	friend std::ostream& operator<<(std::ostream&, const prettyPrint&);
	friend std::istream& operator>>(std::istream&, prettyPrint&);
};

std::ostream& operator<<(std::ostream& os, const prettyPrint& p){
	bool first=true;
	os << '{';
	for(auto x : p.v){
		if(!first)
			os << ", ";
		else
			first=false;
		os << x;
	}
	return(os << '}');
}

std::istream& operator>>(std::istream& is, prettyPrint& p){
	char c;
	is >> c;
	if(c!='{'){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	bool done=false;
	std::vector<double> buffer;
	while(!done){
		c=is.peek();
		if(is.eof())
			break;
		switch(c){
			case '}':
				is >> c;
				done=true;
				break;
			case ',':
				is >> c;
				//fall through
			default:
			{
				double d;
				is >> d;
				if(is.fail())
					return(is);
				buffer.push_back(d);
			}
				break;
		}
	}
	p.v.swap(buffer);
	return(is);
}

///A container to help with parsing index/value pairs for fixing likelihood parameters from the commandline
struct paramFixSpec{
	struct singleParam : std::pair<unsigned int, double>{};
	std::vector<std::pair<unsigned int,double>> params;
};

std::istream& operator>>(std::istream& is, paramFixSpec::singleParam& p){
	char c;
	is >> c;
	if(c!='{'){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	
	is >> p.first >> c >> p.second;
	if(c!=','){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	
	is >> c;
	if(c!='}'){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	return(is);
}
std::istream& operator>>(std::istream& is, paramFixSpec& p){
	char c;
	is >> c;
	if(c!='{'){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	bool done=false;
	while(true){
		c=is.peek();
		if(is.eof())
			break;
		switch(c){
			case '}':
				is >> c;
				done=true;
				break;
			case ',':
				is >> c;
				//fall through
			case '{':
				{
					paramFixSpec::singleParam pa;
					is >> pa;
					if(is.fail())
						return(is);
					p.params.push_back(pa);
				}
				break;
			default:
				is.setstate(is.rdstate()|std::ios::failbit);
				return(is);
		}
	}
	return(is);
}

template<typename HistType>
void accumulateHistogram(HistType& h){
	typename HistType::dataType sum{0};
	for(auto it=h.begin(), end=h.end(); it!=end; it++){
		*it+=sum; //order is important. . .
		sum=*it;
	}
}

template<typename HistType>
void accumulateHistogramReverse(HistType& h){
	typename HistType::dataType sum{0};
	for(auto it=h.rbegin(), end=h.rend(); it!=end; it++){
		*it+=sum; //order is important. . .
		sum=*it;
	}
}

//search for the the value x such that func(c)==target where c is in [a,b]
//and a=func(a), b=func(b), assuming that func is monotonic on this domain.
//This is some sort of simplified form of Brent's method.
template<typename FuncType>
double findRoot(FuncType func, double target, double a, double b, double fa, double fb, double accuracy=1e-6, bool verbose=false){
	bool increasing=fb>fa;
	fa-=target;
	fb-=target;
	
	double c=(a+b)/2, fc;
	double ep=b-a,e=ep,p,q,r,s,t;
	
	while((b-a)>accuracy){
		fc=func(c)-target;
		
		//calculate the quadratic interpolation
		r=fc/fb;
		s=fc/fa;
		t=fa/fb;
		p=s*(t*(r-t)*(b-c)+(r-1)*(c-a));
		q=(r-1)*(s-1)*(t-1);
		
		//collapse the interval
		if((fc>0.0) != increasing){
			a=c;
			fa=fc;
		}
		else{
			fb=fc;
		}
		if(verbose){
			b=c;
			std::cout << " interval now [" << a << ',' << b << ']' << std::endl;
			std::cout << " interpolated guess = " << c+p/q << std::endl;
		}
		
		//accept the interpolation only if it falls within the current boundaries
		if((c+(p/q))>=a && (c+(p/q))<=b){
			//would like to interpolate, but only do so if recent convergence has been suitably rapid
			if(std::abs(p/q) >= 0.5*std::abs(e)){ //it has not; bisect instead
				if(verbose)
					std::cout << " will bisect (too slow)" << std::endl;
				e=ep;
				ep=0.5*(a+b)-c;
				c=0.5*(a+b);
			}
			else{ //it has; use the interpolation
				if(verbose)
					std::cout << " will interpolate" << std::endl;
				c+=p/q;
				e=ep;
				ep=p/q;
			}
		}
		else{ //otherwise, bisect
			if(verbose)
				std::cout << " will bisect (out-of bounds)" << std::endl;
			e=ep;
			ep=0.5*(a+b)-c;
			c=0.5*(a+b);
		}
	}
	return(c);
}

//Given a function which is known to have a root in the domain [a,limit) or (limit,a]
//and an initial guess b which is on the same side of a as the root (>a or <a, respectively)
//finds an interval which strictly brackets the root.
//THe return value is a pair consisting of the other endpoint of the interval, and the function value there
template<typename FuncType>
std::pair<double,double> bracketRoot(FuncType func, double target, double a, double b, double fa, double limit=std::numeric_limits<double>::quiet_NaN()){
	bool dir=b>a;
	if(std::isnan(limit))
		limit=(dir?1:-1)*std::numeric_limits<double>::infinity();
	//std::cout << " bracketing domain is [" << a << ',' << limit << ')' << std::endl;
	fa-=target;
	double c, fc=fa;
	double fact=1;
	while((fc>0)==(fa>0)){
		c=a+fact*(b-a);
		if((dir && c>limit) || (!dir && c<limit)){
			//std::cout << " clipping endpoint to limit" << std::endl;
			c=limit;
		}
		//std::cout << " testing endpoint " << c << std::endl;
		fc=func(c)-target;
		//std::cout << "  (" << fc << ')' << std::endl;
		if(c==limit)
			break;
		fact*=2;
	}
	return(std::make_pair(c,fc+target));
}

///Compute the largest value of the astrophysical flux normalization (assumed to be parameter 2)
///allowed at a given confidence level, estimating significance with Wilks Theorem.
///\param prob the likelihood problem being considered
///\param inputSeed A seed for the unconstrained best fit parameter values. Optional,
///                 but calculation is faster if this has already been computed
///\param indicesToFix The indices of any likelihood parameters which are to be held fixed in
///                    both the constarined and unconstrained fits. Must not contain the index
///                    of the parameter whose limit is being computed.
template<typename ProblemType>
double wilksLim(const ProblemType& prob, double conf,
                boost::optional<std::vector<double>> inputSeed=boost::optional<std::vector<double>>(),
                std::vector<unsigned int> indicesToFix={}){
	std::vector<double> seed=inputSeed?*inputSeed:prob.getSeed();
	fitResult best=doFitLBFGSB(prob,seed,indicesToFix);
	
	std::cout << " best fit is: ";
	for(unsigned int i=0; i<best.params.size(); i++)
		std::cout << best.params[i] << ' ';
	std::cout << std::setprecision(10) << best.likelihood << std::setprecision(6) << ' ' << (best.succeeded?"succeeded":"failed") << std::endl;
	
	double targetR=boost::math::quantile(boost::math::chi_squared_distribution<double>(1),conf);
	std::cout << " quantile for a confidence level of " << conf << " is " << targetR << std::endl;
	
	indicesToFix.push_back(2);
	//exponential bracketing
	std::cout << " bracketing. . . " << std::endl;
	double sigMax=std::max(best.params[2],1.0);
	double sigMin=best.params[2];
	double rSigMin=0, rSigMax;
	double R=0;
	while(R<targetR){
		sigMax*=2;
		std::cout << "  testing " << sigMax << std::endl;
		//evaulate R
		std::vector<double> point=best.params;
		point.resize(6);
		point[2]=sigMax;
		fitResult result=doFitLBFGSB(prob,point,indicesToFix);
		
		std::cout << "  fit is ";
		for(unsigned int i=0; i<result.params.size(); i++)
			std::cout << result.params[i] << ' ';
		std::cout << std::setprecision(10) << result.likelihood << std::setprecision(6) << ' ' << (result.succeeded?"succeeded":"failed") << std::endl;
		
		R=-2*(best.likelihood-result.likelihood); //likelihoods are already logarithms
		std::cout << "  R=" << R << std::endl;
		if(R<targetR){
			sigMin=sigMax;
			rSigMin=R;
		}
	}
	rSigMax=R;
	
	std::cout << " locating root. . . " << std::endl;
	std::vector<double> point=best.params;
	auto bestLikelihood=best.likelihood;
	auto testFunc=[&prob,&indicesToFix,point,bestLikelihood](double test) mutable->double{
		std::cout << "  testing " << test << std::endl;
		point[2]=test;
		fitResult result=doFitLBFGSB(prob,point,indicesToFix);
		
		std::cout << "  fit is ";
		for(unsigned int i=0; i<result.params.size(); i++)
			std::cout << result.params[i] << ' ';
		std::cout << std::setprecision(10) << result.likelihood << std::setprecision(6) << ' ' << (result.succeeded?"succeeded":"failed") << std::endl;
		
		double R=-2*(bestLikelihood-result.likelihood); //likelihoods are already logarithms
		std::cout << "  R=" << R << std::endl;
		return(R);
	};
	double limit=findRoot(testFunc, targetR, sigMin, sigMax, rSigMin, rSigMax, /*8e-3*/1e-5);
	std::cout << " Upper limit: " << limit << std::endl;
	return(limit);
}

///Compute (double) the likelihood ratio between the global best fit with all parameters free,
///to the constarined fit with some parameters fixed to zero
///\param prob The likelihood problem for which to evaluate the ratio
///\param paramsToFix The parameters to fix and the values to which to fix them for the constrained fit
///\param indicesToFix The indices of parameters which should be held fixed in both fits
///\param inputSeed A set of seed values for the parameter which is or is close to the global minimum
template<typename ProblemType>
double computeLikelihoodRatio(const ProblemType& prob, std::vector<std::pair<unsigned int,double>> paramsToFix={},
							  std::vector<unsigned int> indicesToFix={},
							  boost::optional<std::vector<double>> inputSeed=boost::optional<std::vector<double>>()){
	std::vector<double> seed=inputSeed?*inputSeed:prob.getSeed();
	fitResult best=doFitLBFGSB(prob,seed,indicesToFix);
	
	if(!quiet){
		std::cout << " best fit is: ";
		for(unsigned int i=0; i<best.params.size(); i++)
			std::cout << best.params[i] << ' ';
		std::cout << std::setprecision(10) << best.likelihood << std::setprecision(6) << ' ' << (best.succeeded?"succeeded":"failed") << std::endl;
	}
	
	for(auto fixParam : paramsToFix){
		seed[fixParam.first]=fixParam.second;
		indicesToFix.push_back(fixParam.first);
	}
	fitResult constrained=doFitLBFGSB(prob,seed,indicesToFix);
	
	if(!quiet){
		std::cout << " constrained fit is: ";
		for(unsigned int i=0; i<constrained.params.size(); i++)
			std::cout << constrained.params[i] << ' ';
		std::cout << std::setprecision(10) << constrained.likelihood << std::setprecision(6) << ' ' << (constrained.succeeded?"succeeded":"failed") << std::endl;
	}
	
	double R=-2*(best.likelihood-constrained.likelihood); //likelihoods are already logarithms
	if(R<0)
		R=0;
	return(R);
}

///Estimate the significance of the best fit astrophysical flux compared to an atmospheric
///only hypothesis using Wilks Theorem
///\param prob The likelihood problem for which to estimate the significance
///\param inputSeed A set of seed values for the parameter which is or is close to the
///                 unconstrained best fit. This is optional but makes the calculation
///                 more efficient.
///\param indicesToFix The indices of parameters which should be held fixed in both fits
template<typename ProblemType>
double wilksNonZeroSig(const ProblemType& prob,
                       boost::optional<std::vector<double>> inputSeed=boost::optional<std::vector<double>>(),
                       std::vector<unsigned int> indicesToFix={}){
	double R=computeLikelihoodRatio(prob,{{2,0.0}},indicesToFix,inputSeed);
	
	if(R<0){ //sometimes the minimizer stops just a little too early, yielding slightly corrupt likelihood ratios
		R=0; //just paper over it, the difference isn't particularly important in practice
		std::cout << "P=" << 1 << " (" << -10 << " sigma)" << std::endl;
		return(-10);
	}
	std::cout << " R=" << R << std::endl;
	auto dist=boost::math::chi_squared_distribution<double>(1);
	double intProb=boost::math::cdf(dist,R);
	boost::math::normal_distribution<double> gauss;
	if(intProb>=1.0){
		std::cout << "probability overflow, rounding to 10 sigma" << std::endl;
		return(10.);
	}
	//double significance=boost::math::quantile(gauss,(1+intProb)/2); //two-tailed
	double significance=boost::math::quantile(gauss,intProb); //one-tailed
	std::cout << "P=" << 1-intProb << " (" << significance << " sigma)" << std::endl;
	return(significance);
}

///Sample test statistics (likelihood ratios) to compute the exact significance of the best
///fit astrophysical flux compared to an atmospheric only hypothesis
///\param simulation1 A set of simulated events which will be used for fitting
///\param simulation2 A set of simulated events from which trial realizations will be drawn
///\param histTemplate A histogram which will be used as the template for the binning in the
///                    likelihood fit
///\param binner An object which will put events into the fit histograms
///\param priors The set of priors to use for the likelihood
///\param wm The object which will compute the weights of simulated events for a given hypothesis
///\param rng The random engine for sampling realizations
///\param testPoint The point in the model space at which the realizations will be sampled
///\param expectedEvents The number of events to sample in each realization
///\param nTests The number of realizations to sample and test
///\param fixedParams The parameters
template<typename HistType, typename BinnerType, typename PriorsType, typename WeighterMaker, typename RNG>
void FCNonZeroSig(const std::deque<Event>& simulation1, const std::deque<Event>& simulation2,
                  HistType histTemplate, BinnerType binner,
                  PriorsType priors, WeighterMaker wm,
                  std::vector<double> seed, RNG& rng,
				  std::vector<double> testPoint, unsigned int expectedEvents,
				  unsigned int nTests=10, paramFixSpec fixedParams=paramFixSpec()){
	std::vector<std::deque<Event>> allSimSets{simulation1};
	std::vector<HistType> allSimHists;
	allSimHists.push_back(makeEmptyHistogramCopy(histTemplate));
	bin(simulation1,allSimHists.back(),binner);
	
	auto weighter=wm(testPoint);
	
	std::vector<unsigned int> fixedIndices;
	for(const auto pf : fixedParams.params){
		seed[pf.first]=pf.second;
		fixedIndices.push_back(pf.first);
	}
	
	std::vector<double> testStatistics;
	for(unsigned int i=0; i<nTests; i++){
		auto sample=likelihood::generateSample(weighter,simulation2,expectedEvents,rng);
		if(!quiet)
			std::cout << " sampled " << sample.size() << " events" << std::endl;
		HistType sampleHist=makeEmptyHistogramCopy(histTemplate);
		bin(sample,sampleHist,binner);
		
		using namespace likelihood;
		auto prob=makeLikelihoodProblem<std::reference_wrapper<const Event>,3,8>(sampleHist, allSimHists, priors, {1.0}, simpleDataWeighter(), wm, poissonLikelihood(), seed);
		prob.setEvaluationThreadCount(evalThreads);
		
		double R=computeLikelihoodRatio(prob,{{2,0.0}},fixedIndices);
		std::cout << " R=" << R << std::endl;
		testStatistics.push_back(R);
	}
}

//compute the median upper limit at confidence `confidence` of `nTests` MC trials
template<typename HistType, typename BinnerType, typename PriorsType, typename WeighterMaker, typename RNG>
double computeSensitivity(const std::deque<Event>& simulation1, const std::deque<Event>& simulation2, HistType histTemplate, BinnerType binner, PriorsType priors, WeighterMaker wm, const std::vector<double>& seed, RNG& rng, unsigned int nTests=100, double confidence=0.9){
	std::vector<std::deque<Event>> allSimSets{simulation1};
	std::vector<HistType> allSimHists;
	allSimHists.push_back(makeEmptyHistogramCopy(histTemplate));
	bin(simulation1,allSimHists.back(),binner);
	
	auto weighter=wm(std::vector<double>{0.943,0,0,-0.027,0.186,1.146});
	double expected=0;
	auto testHist=makeEmptyHistogramCopy(histTemplate);
	bin(simulation2,testHist,binner);
	for(typename HistType::dataType bin : testHist){
		for(const Event& event : bin){
			auto w=weighter(event);
			expected+=w;
		}
	}
	
	std::vector<double> limits;
	for(unsigned int i=0; i<nTests; i++){
		auto sample=likelihood::generateSample(weighter,simulation2,expected,rng);
		if(!quiet)
			std::cout << " sampled " << sample.size() << " events" << std::endl;
		HistType sampleHist=makeEmptyHistogramCopy(histTemplate);
		bin(sample,sampleHist,binner);
		
		using namespace likelihood;
		auto prob=makeLikelihoodProblem<std::reference_wrapper<const Event>,3,8>(sampleHist, allSimHists, priors, {1.0}, simpleDataWeighter(), wm, poissonLikelihood(), seed);
		prob.setEvaluationThreadCount(evalThreads);
		
		limits.push_back(wilksLim(prob,confidence));
	}
	
	if(nTests%2==0)
		return((limits[nTests/2]+limits[nTests/2+1])/2);
	return(limits[nTests/2]);
}

template<typename HistType, typename BinnerType, typename PriorsType, typename WeighterMaker, typename RNG>
double computeSensitiveRange(const std::deque<Event>& simulation1, const std::deque<Event>& simulation2, HistType histTemplate, BinnerType binner, PriorsType priors, WeighterMaker wm, RNG& rng, double initialValue, bool upperLimit){
	auto sensFunc=[&](double limit)->double{
		if(upperLimit)
			histTemplate.getAxis(0)->setUpperLimit(limit);
		else
			histTemplate.getAxis(0)->setLowerLimit(limit);
		double sens=computeSensitivity(simulation1,simulation2,histTemplate,binner,priors,wm,rng,100);
		std::cout << "senstivity for " << limit << " is " << sens << std::endl;
		return(sens);
	};
	double baseSens=sensFunc(initialValue);
	double targetSens=1.05*baseSens, boundSens=baseSens;
	
	double bound=initialValue;
	while(boundSens<targetSens){
		if(upperLimit)
			bound/=2;
		else
			bound*=2;
		boundSens=sensFunc(bound);
	}
	if(bound>initialValue){
		std::swap(bound,initialValue);
		std::swap(boundSens,baseSens);
	}
	
	return(findRoot(sensFunc, targetSens, bound, initialValue, boundSens, baseSens, 1e-2));
}

template<typename ProblemType>
double saturatedPoisson(const ProblemType& prob, const std::vector<double>& params){
	auto spProb=prob.makeAlternateLikelihood(likelihood::saturatedPoissonLikelihood());
	return(-spProb.evaluateLikelihood(params));
}

struct timer{
private:
	std::chrono::high_resolution_clock::time_point t1, t2;
	bool running=false;
	std::string name;
	friend std::ostream& operator<<(std::ostream& os, timer& t);
public:
	timer():running(false){}
	explicit timer(std::string name):running(false),name(name){}
	
	~timer(){
		if(running)
			std::cout << *this << std::endl;
	}
	
	std::string getName() const{
		return(name);
	}
	void setName(std::string name){
		this->name=name;
	}
	
	void start(){
		if(!running){
			t1 = std::chrono::high_resolution_clock::now();
			running=true;
		}
	}
	void start(std::string name){
		if(!running){
			this->name=name;
			t1 = std::chrono::high_resolution_clock::now();
			running=true;
		}
	}
	void stop(){
		if(running){
			t2 = std::chrono::high_resolution_clock::now();
			running=false;
		}
	}
	double get(){
		if(running)
			stop();
		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		return(time_span.count());
	}
};

std::ostream& operator<<(std::ostream& os, timer& t){
	if(quiet)
		return(os);
	if(!t.getName().empty())
		os << t.getName() << ": ";
	return(os << t.get() << " seconds");
}

//precompute weights as much as possible
/////////////////////////////////////////////////// VERY IMPORTANT // MUY IMPORTANTE ///////////////////////////////////////////////////
template<typename ContainerType, typename WeighterType>
void initializeSimulationWeights(ContainerType& simulation, const WeighterType& convPionWeighter, const WeighterType& convKaonWeighter, const WeighterType& promptWeighter, const OversizeWeighter& osw){
	using iterator=typename ContainerType::iterator;
	auto cache=[&](iterator it, iterator end){
		for(; it!=end; it++){
			auto& e=*it;
      LW::Event lw_e {e.leptonEnergyFraction,
                e.injectedEnergy,
                e.totalColumnDepth,
                e.inelasticityProbability,
                e.intX,
                e.intY,
                e.injectedMuonEnergy,
                e.injectedMuonZenith,
                e.injectedMuonAzimuth,
                static_cast<particleType>(e.primaryType),
                e.year};

			double osweight = osw.EvaluateOversizeCorrection(e.energy, e.zenith);
      e.cachedConvPionWeight=convPionWeighter(lw_e)*e.cachedLivetime*osweight;
			e.cachedConvKaonWeight=convKaonWeighter(lw_e)*e.cachedLivetime*osweight;
			//e.cachedPromptWeight=promptWeighter(lw_e)*e.cachedLivetime;
			e.cachedPromptWeight=0.;
     // std::cout << osweight << " " << e.cachedConvPionWeight << " " << e.cachedConvKaonWeight << " " << e.cachedPromptWeight << std::endl;
      // we will set this to zero. Love. CA.
			e.cachedAstroWeight=0.;
		}
	};

	ThreadPool pool(evalThreads);
	iterator it=simulation.begin(), end=simulation.end();
	while(true){
		unsigned long dist=std::distance(it,end);
		if(dist>=1000UL){
			pool.enqueue(cache,it,it+1000);
			it+=1000;
		}
		else{
			pool.enqueue(cache,it,end);
			break;
		}
	}
}


int main(int argc, char* argv[]){
  std::string compact_data_path =        "../compact_data/";
  std::string plot_path =                "../plots/";
  std::string output_path =              "../output/";
  //std::string squids_files_path =        "/data/ana/NuFSGenMC/IC86_hypotheses/fluxes/flux_ln0/";
  std::string squids_files_path =        "/data/ana/NuFSGenMC/IC86_hypotheses/fluxes_new/flux_ln0/";
  std::string prompt_squids_files_path = "/data/ana/NuFSGenMC/IC86_hypotheses/fluxes/prompt_0/";
  std::string xs_spline_path =           "/data/ana/NuFSGenMC/CrossSections/";
  std::string data_path =                "/data/ana/NuFSGenMC/Data/";
  std::string mc_path =                  "/data/ana/NuFSGenMC/Merged/";
  //std::string oversize_function_main =   "/data/ana/NuFSGenMC/OversizeCorrections/NullCorrection.dat";
  //std::string oversize_function_dc =     "/data/ana/NuFSGenMC/OversizeCorrections/NullCorrection.dat";
  std::string oversize_function_main =   "/data/ana/NuFSGenMC/OversizeCorrections/os_corr_0.99.dat";
  std::string oversize_function_dc =     "/data/ana/NuFSGenMC/OversizeCorrections/os_corr_0.99.dat";
  std::string domeff_spline_path =       "/data/ana/NuFSGenMC/DomEffSplines/";

  std::string model_name = "";
  std::string delta_model_name = "Honda+Gaisser";
	uint32_t rngSeed=0;
	double minFitEnergy=4e2;
	//double minFitEnergy=7e2;
	//double maxFitEnergy=8.e3;
	double maxFitEnergy=2.e4;
  double minCosth = -1.;
  double maxCosth = 0.2;
	double minAstroEnergy=0.0;
	double maxAstroEnergy=1e10;
  double minAzimuth = 0.;
  double maxAzimuth = 2.*pi<double>();
  double minCOGZ = std::numeric_limits<double>::min();
  double maxCOGZ = std::numeric_limits<double>::max();
	bool doComputeSensitivity=false;
	unsigned int sensitivityTrials=100;
	int sampleTestStatistics=0;
	std::string outputFile;
	bool drawPlots=true;
	paramFixSpec fixedParams;
	std::vector<double> existingBest;
	std::vector<double> fitSeed{1.02,0,0,0.05,.0985,1.1,1,0};
  std::vector<double> dataChallenge_nuisance_parameters{1,0,0,0,.1,1,1,0};
  std::vector<double> data{1,0,0,0,.1,1,1,0};
	prettyPrint existingBest_(existingBest);
	prettyPrint fitSeed_(fitSeed);
	prettyPrint dataChallenge_nuisance_parameters_(dataChallenge_nuisance_parameters);
  unsigned int number_of_data_challenges = 1;
  std::vector<double> dataChallenge_seeds {1234};
  prettyPrint dataChallenge_seeds_(dataChallenge_seeds);
  bool use_datachallenge_histogram=false;
	bool writeCompact=false;
	bool readCompact=false;
	bool computeWilksSig=false;
	bool doDataChallenge=false;
	bool exitAfterLoading=false;
	int yearsIC86=1;
  bool UseBurnsample=true;
	bool computeParameterErrors=false;
  double th24 = 0;
  double dm41sq = 0;
  double th24_null = 0;
  double dm41sq_null = 0;
  bool use_gsl_minimizer = false;
  bool dump_data = false;
  bool dump_mc_data = false;
  bool dump_real_data = false;
  bool dump_fit = false;
  bool save_flux = false;
  bool save_dc_sample = false;
  bool doBootstrap = false;
  std::string simulation_to_load = "nufsgen_mie_0_99";
  std::string dcsimulation_to_load = "nufsgen_mie_0_99";
  std::string xs_model_name = "";
  std::string data_challenge_histogram_filepath = "";

	OptionParser op;
	op.addOption({"j","numEvalThreads"},evalThreads,"Number of threads to use when computing likelihoods. Must be >0.");
	op.addOption({"s","rngSeed"},rngSeed,"Seed value for the RNG, if 0 system clock time is used.");
	op.addOption({"q","quiet"},[&](){quiet=true;},"Reduce chatter to stdout.");
	op.addOption("minFitEnergy",minFitEnergy,"Minimum reconstructed energy in GeV for events to participate in the likelihood fit.");
	op.addOption("maxFitEnergy",maxFitEnergy,"Maximum reconstructed energy in GeV for events to participate in the likelihood fit.");
	op.addOption("minCosth",minCosth,"Minimum reconstructed cos(th) for events to participate in the likelihood fit.");
	op.addOption("maxCosth",maxCosth,"Maximum reconstructed cos(th) for events to participate in the likelihood fit.");
	op.addOption("minAstroEnergy",minAstroEnergy,"Minimum true energy in GeV for the astrophysical flux to be non-zero.");
	op.addOption("maxAstroEnergy",maxAstroEnergy,"Maximum true energy in GeV for the astrophysical flux to be non-zero.");
	op.addOption("computeSensitivity",doComputeSensitivity,"Whether to compute a sensitivity estimate using MC trials.");
	op.addOption("sensitivityTrials",sensitivityTrials,"Number of trials to use when computing a sensitivity estimate.");
	op.addOption("sampleTestStatistics",sampleTestStatistics,"Whether to sample likelihood ratio test statistics using MC trials.");
	op.addOption({"o","output"},outputFile,"If specified as non-empty the file into which to direct standard output.");
	op.addOption({"p","drawPlots"},drawPlots,"Whether to generate plots directly.");
	op.addOption({"f","fix"},fixedParams,"Parameters to hold fixed to particular values in the fit.");
	op.addOption("bestFit",existingBest_,"Specifiy a best fit computed in a previous run so that it need not be recomputed.");
	op.addOption("fitSeed",fitSeed_,"Specifiy the set of parameters to use when seeding the fitter.");
	op.addOption("dataChallenge_nuisance_parameters",dataChallenge_nuisance_parameters_,"Specifiy the set of nuisance parameters to use when performing a data challenge.");
	op.addOption("writeCompact",[&](){writeCompact=true;},"Write all input data back out in the compact format.");
	op.addOption("readCompact",[&](){readCompact=true;},"Read all input data in the compact format\n    "
				 "(requires that the same version of the program has been run previously with --writeCompact).");
	op.addOption("exitAfterLoading",[&](){exitAfterLoading=true;},"Exit immediately after writing compact data.");
  op.addOption("UseBurnsample",[&](){UseBurnsample = true;},"If True burn sample will be use, else full sample.");
  op.addOption("UseFullsample",[&](){UseBurnsample = false;},"If True full sample will be use.");
  op.addOption("DumpDCFit",[&](){dump_fit = true;},"If True data challenge fit for the H1 hypohteiss will be dump.");
  op.addOption("DumpDCData",[&](){dump_data = true;},"If True data challenge realization will be dump.");
  op.addOption("DumpData",[&](){dump_real_data = true;},"If True data realization will be dump.");
  op.addOption("DumpMC",[&](){dump_mc_data = true;},"If True the MC will be dump.");
	op.addOption("wilksSig",[&](){computeWilksSig=true;},"Compute a significance for a nonzero astrophysical flux with Wilks theorem.");
	op.addOption("dataChallenge",[&](){doDataChallenge=true;},"Run fit trials on MC realizations.");
	op.addOption("save_flux",[&](){save_flux=true;},"Save flux..");
	op.addOption("yearsIC86",yearsIC86,"Number of years of IC86 equivalent livetime.");
	op.addOption("parameterErrors",[&](){computeParameterErrors=true;},"Compute 1 sigma Wilks Theorem error estimates on all free fit parameters.");
  op.addOption("th24_null",th24_null,"th24 sterile neutrino mixing paramater null hypothesis [rad].");
  op.addOption("th24",th24,"th24 sterile neutrino mixing paramater [rad].");
  op.addOption("dm41sq_null",dm41sq_null,"dm41sq sterile neutrino mass paramater null hypothesis[ev^2].");
  op.addOption("dm41sq",dm41sq,"dm41sq sterile neutrino mass paramater [ev^2].");
  op.addOption("simulation_to_load",simulation_to_load, "String that specifies the simulation to load");
  op.addOption("dataChallenge_simulation_to_load",dcsimulation_to_load, "String that specifies the data challenge null hypothesis simulation to load");
  op.addOption("plot_path",plot_path, "String that specifies the plot data path");
  op.addOption("compact_data_path",compact_data_path, "String that specifies the compact data path");
  op.addOption("data_path",data_path, "String that specifies the data path");
  op.addOption("mc_path",mc_path, "String that specifies the mc path");
  op.addOption("squids_files_path",squids_files_path, "String that specifies the path to the squids files");
  op.addOption("prompt_squids_files_path",prompt_squids_files_path, "String that specifies the path to the prompt squids files");
  op.addOption("xs_spline_path",xs_spline_path, "String that specifies the path to the cross section splines");
  op.addOption("xs_model_name",xs_model_name,"Name of the cross section model to use");
  op.addOption("domeff_spline_path",domeff_spline_path,"Path to the DOM efficiencies splines");
  op.addOption("model_name",model_name, "Hadronic model name, if empty will use Honda+Gaisser");
  op.addOption("delta_model_name",delta_model_name, "Hadronic model name, if empty will use Honda+Gaisser");
  op.addOption("output_path",output_path, "String that specifies the output data path");
  op.addOption("use_gsl_minimizer",[&](){use_gsl_minimizer = true;},"If use_gsl_minimizer is set to true then GSL minimizer will be used.");
  op.addOption("doBootstrap",[&](){doBootstrap = true;},"If doBootstrap then bootstrap will be performed.");
  op.addOption("oversize_function_main",oversize_function_main,"Main oversize correction function to use");
  op.addOption("oversize_function_dc",oversize_function_dc,"Data challenge oversize correction function to use");
  op.addOption("number_of_data_challenges",number_of_data_challenges,"Integer specifying the number of data challenges to perform.");
  op.addOption("dataChallenge_seeds",dataChallenge_seeds_,"Seeds used for each the data challenges.");
  op.addOption("minAzimuth",minAzimuth,"Minimul azimuth value [radians].");
  op.addOption("maxAzimuth",maxAzimuth,"Minimul azimuth value [radians].");
  op.addOption("minCOGZ",minCOGZ,"Minimul min COGZ value [meters].");
  op.addOption("maxCOGZ",maxCOGZ,"Minimul max COGZ value [meters].");
  op.addOption("use_datachallenge_histogram",[&](){use_datachallenge_histogram = true;},"If use_datachallenge_histogram is set to true then the DC MC would not be loaded and the histogram will be loaded instead.");
  op.addOption("data_challenge_histogram_filepath",data_challenge_histogram_filepath,"Path to the histogram use for data challenge when use_datachallenge_histogram is set to true.");
  op.addOption("save_datachallenge_histogram",[&](){save_dc_sample = true;},"If save_dc_sample is set to true then the DC MC will be saved with filename data_challenge_histogram_filepath.");

	auto args=op.parseArgs(argc,argv);
	if(op.didPrintUsage())
		return(0);
	
	if(evalThreads==0){
		std::cerr << "Requesting likelihood evaluation with 0 threads of execution makes no sense" << std::endl;
		return(1);
	}

  if((number_of_data_challenges != dataChallenge_seeds.size()) and doDataChallenge){
		std::cerr << "Number of data challenges does not match data challenges seeds." << std::endl;
		return(1);
  }

	if(exitAfterLoading && !writeCompact){
		std::cerr << "--exitAfterLoading does not make sense without --writeCompact" << std::endl;
		return(1);
	}
	if(!existingBest.empty() && existingBest.size()!=7){
		std::cerr << "Existing fit (--bestFit) must have 7 parameter values" << std::endl;
		return(1);
	}
	if(fitSeed.size()!=8){
		std::cerr << "Fit seed must have exactly 8 parameter values" << std::endl;
		std::cerr << fitSeed << std::endl;
		return(1);
	}
	
	std::ofstream altOutput;
	std::streambuf* stdoutBackup;
	if(!outputFile.empty()){
		altOutput.open(outputFile.c_str());
		if(!altOutput.good()){
			std::cerr << "Failed to open log file " << outputFile << " for writing" << std::endl;
			return(1);
		}
		stdoutBackup=std::cout.rdbuf();
		std::cout.rdbuf(altOutput.rdbuf());
	}

  // read dom efficincies splines
  if(!quiet)
    std::cout << "Begin Spline Loading" << std::endl;
  domEffConv2011 = std::unique_ptr<Splinetable>(new Splinetable(domeff_spline_path+"/conv_IC86.fits"));
  domEffPrompt2011 = std::unique_ptr<Splinetable>(new Splinetable(domeff_spline_path+"/prompt_IC86.fits"));
  //domEffConvPion2011 = std::unique_ptr<Splinetable>(new Splinetable(domeff_spline_path+"/conv_pion_IC86.fits"));
  //domEffConvKaon2011 = std::unique_ptr<Splinetable>(new Splinetable(domeff_spline_path+"/conv_kaon_IC86.fits"));
  if(!quiet)
    std::cout << "End Spline Loading" << std::endl;
	
	OversizeWeighter osw_main( oversize_function_main ) ;

	const std::string mc_dataPath= mc_path;
	const std::string data_dataPath= data_path;
	std::deque<Event> mainSimulation;
	std::deque<Event> burnsample;
	
	//livetime for each year of data taking
  std::map<int,double> livetime;
  if (UseBurnsample)
    livetime = std::map<int,double> {{2011,758.59*60*60}}; // IC86 burnsample only
  else
    livetime = std::map<int,double> {{2011,8249.6*3600}}; // IC86 burnsample only

  // number of events per file
  const std::map<std::string,run> simInfo = {
                                {"nufsgen_mie_0_9",run(mc_dataPath,"noover_mie_0.9.h5",particleType::NuMu,0.9,
                                        LW::SimDetails {1,615360000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_mie_0_95",run(mc_dataPath,"noover_mie_0.95.h5",particleType::NuMu,0.95,
                                        LW::SimDetails {1,292140000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_mie_0_99",run(mc_dataPath,"noover_mie_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,734400000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_mie_1_089",run(mc_dataPath,"noover_mie_1.089.h5",particleType::NuMu,1.089,
                                        LW::SimDetails {1,584640000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_mie_1_1979",run(mc_dataPath,"noover_mie_1.1979.h5",particleType::NuMu,1.1979,
                                        LW::SimDetails {1,440400000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_mie_0_99_sel",run(mc_dataPath,"noover_mie_0.99_sel.h5",particleType::NuMu,0.99,
                                                        LW::SimDetails {1,57300000, //total events                                            
                                                        2011,//year                                                                              
                                                        800, //injection radius                                                                   
                                                        0,2*pi<double>(), //azimuth range                                                         
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range                                      
                                                        2e2,1e6, //energy range                                                                   
                                                        2} // spectral index                                                                      
                                )},
                                {"nufsgen_lea_0_9",run(mc_dataPath,"noover_lea_0.9.h5",particleType::NuMu,0.9,
                                        LW::SimDetails {1,39416000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                /*{"nufsgen_lea_0_99",run(mc_dataPath,"noover_lea_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,39568000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},*/
                                {"nufsgen_lea_1_1979",run(mc_dataPath,"noover_lea_1.1979.h5",particleType::NuMu,1.1979,
                                        LW::SimDetails {1,39616000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_lea_0_99",run(mc_dataPath,"lea_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,884940000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_noholeice_mie_0_99",run(mc_dataPath,"noholeice_mie_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,295710000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_ice1_mie_0_99",run(mc_dataPath,"ice1_mie_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,298260000,//total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_ice2_mie_0_99",run(mc_dataPath,"ice2_mie_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,298440000,//total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )}
                              };

	timer t("Data loading");
	t.start();

  const std::vector<std::string> simSetsToLoad = { simulation_to_load };

	if(readCompact){
		try{
			unsplatData(compact_data_path+"/"+simulation_to_load+"_compact_data.dat",getFileChecksum(argv[0]),burnsample,mainSimulation);
			if(!quiet){
				std::cout << "Loaded " << burnsample.size() << " experimental events" << std::endl;
				std::cout << "Loaded " << mainSimulation.size() << " events in main simulation set" << std::endl;
			}
		} catch(std::runtime_error& re){
			std::cerr << re.what() << std::endl;
			std::cerr << "Failed to load compact data" << std::endl;
			return(1);
		}
	}
	else{
		try{
			burnsample=loadExperimentalData(data_dataPath,UseBurnsample);
		} catch(std::exception& ex){
			std::cerr << "Problem loading experimental data: " << ex.what() << std::endl;
			return(1);
		}
		if(!quiet)
			std::cout << "Loaded " << burnsample.size() << " experimental events" << std::endl;

		const bool loadTargeted=true;
		try{
      loadSimulatedData(mainSimulation,mc_dataPath,livetime,simInfo,simSetsToLoad,loadTargeted);
		} catch(std::exception& ex){
			std::cerr << "Problem loading simulated data: " << ex.what() << std::endl;
			return(1);
		}
		if(!quiet)
			std::cout << "Loaded " << mainSimulation.size() << " events in main simulation set" << std::endl;
	}

	if(writeCompact){
		try{
			splatData(compact_data_path+"/"+simulation_to_load+"_compact_data.dat",
                getFileChecksum(argv[0]),burnsample,mainSimulation);
		} catch(std::runtime_error& re){
			std::cerr << re.what() << std::endl;
			std::cerr << "Failed to save compact data" << std::endl;
			return(1);
		}
		if(exitAfterLoading)
			return(0);
	}

  t.stop();
  if(!quiet)
    std::cout << t << std::endl;

  if(dump_real_data){
    std::ofstream dump_data;
    if(UseBurnsample)
     dump_data = std::ofstream(output_path + "/" +"IC86burnsample.dat");
    else
     dump_data = std::ofstream(output_path + "/" +"IC86fullsample.dat");

    dump_data << std::scientific;
    dump_data.precision(6);
    dump_data.fill('0');
    for(Event& e : burnsample) {
      //if( e.energy > 4.e2 and e.energy < 2.e4 ){
      //if( e.energy > 4.e2 and e.energy < 2.e4 ){
      if( true ){
      //{
        //dump_data << cos(e.zenith) << ' ' << e.azimuth << ' ' << e.energy << ' ';
        dump_data << e.runid << ' ' << e.eventid << ' ';
        dump_data << e.energy << ' ' << cos(e.zenith) << ' ';
#ifdef _FULL_EVENT_
        dump_data << e.cogx << ' ' << e.cogy << ' ' << e.cogz << ' ';
        dump_data << e.rlogl << ' ';
        dump_data << e.time << ' ';
        dump_data << e.n_chan << ' ' << e.n_chan_nodc << ' ';
#endif
        //dump_data << 1. << std::endl;
        dump_data << std::endl;
      }
    }

    dump_data.close();
  }
  //exit(0);

#ifdef _FULL_EVENT_
  // hack to do split fits
  for(auto it = mainSimulation.begin(); it!=mainSimulation.end(); it++){
    if ((*it).azimuth < minAzimuth or (*it).azimuth > maxAzimuth)
      mainSimulation.erase(it);
    else if ((*it).cogz < minCOGZ or (*it).cogz > maxCOGZ )
      mainSimulation.erase(it);
  }

  for(auto it = burnsample.begin(); it!=burnsample.end(); it++){
    if ((*it).azimuth < minAzimuth or (*it).azimuth > maxAzimuth)
      burnsample.erase(it);
    else if ((*it).cogz < minCOGZ or (*it).cogz > maxCOGZ )
      burnsample.erase(it);
  // lala
  if(!quiet){
  }
    std::cout << "Performed COGZ/Azimuth cuts." << std::endl;
    std::cout << "Loaded " << mainSimulation.size() << " events in main simulation set" << std::endl;
    std::cout << "Loaded " << burnsample.size() << " events in data set" << std::endl;
  }
#endif

/////////////////////////////////////////////////// VERY IMPORTANT // MUY IMPORTANTE ///////////////////////////////////////////////////
// HERE WE CONSTRUCT THE LW WEIGHTERS FOR THE OSCILLATION HYPOTHESIS AND WEIGH THE MC //
// AQUI CONSTRUIMOS LOS PESADORES DE LW PARA UNA HIPOTESIS DE OSCILACION Y PESAMOS EL MC //
	t.start("WeightingMC");
  // flux weighters // pesadores del flujo
  if(!quiet)
    std::cout << "Begin Loading nuSQuIDS objects." << std::endl;
  std::shared_ptr<LW::SQUIDSFlux> flux_kaon,flux_pion;
  if (model_name == ""){
    flux_kaon = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "kaon_atmospheric_"+std::to_string(dm41sq)+"_"+std::to_string(th24)+".hdf5");
    flux_pion = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "pion_atmospheric_"+std::to_string(dm41sq)+"_"+std::to_string(th24)+".hdf5");
  } else {
    flux_kaon = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "kaon_atmospheric_"+std::to_string(dm41sq)+"_"+std::to_string(th24)+"_"+model_name+".hdf5");
    flux_pion = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "pion_atmospheric_"+std::to_string(dm41sq)+"_"+std::to_string(th24)+"_"+model_name+".hdf5");
  }

  // prompt things
  //std::shared_ptr<LW::SQUIDSFlux> flux_prompt = std::make_shared<LW::SQUIDSFlux>(prompt_squids_files_path + "prompt_atmospheric_"+std::to_string(dm41sq)+"_"+std::to_string(th24)+".hdf5");
  // assume that prompt is O(NH) prompt
  std::shared_ptr<LW::SQUIDSFlux> flux_prompt = std::make_shared<LW::SQUIDSFlux>(prompt_squids_files_path + "prompt_atmospheric_0.000000_0.000000.hdf5");

  // cross section weighter // pesador de seccion de choque
  std::shared_ptr<LW::CrossSectionFromSpline> xsw = std::make_shared<LW::CrossSectionFromSpline>(static_cast<std::string>(xs_spline_path),xs_model_name);
  // mc generation weighter // pesador de la generacion del mc
  LW::mcgenWeighter mcw;
  for( std::string sim_name : simSetsToLoad )
    mcw.addGenerationSpectrum(simInfo.find(sim_name)->second.details);
  // lw // lw
  LW::LeptonWeighter PionFluxWeighter(flux_pion,xsw,mcw);
  LW::LeptonWeighter KaonFluxWeighter(flux_kaon,xsw,mcw);
  LW::LeptonWeighter PromptFluxWeighter(flux_prompt,xsw,mcw);
  // initialize simulation weights // inicializa los pesos de la simulacion
	initializeSimulationWeights(mainSimulation,PionFluxWeighter,KaonFluxWeighter,PromptFluxWeighter,osw_main);

  if(save_flux)
  {
    std::cout << "Saving flux model" << std::endl;
    std::string sterile_params_str = "dm41sq_"+ std::to_string(dm41sq)+"_th24_"+std::to_string(th24);
    std::ofstream flux_file(output_path + simulation_to_load + "_flux_sampled_"+
                                  sterile_params_str+"_"+model_name+"_.dat");
    for(auto& e : mainSimulation){
        flux_file << static_cast<int>(e.primaryType) << ' ';
        flux_file << cos(e.injectedMuonZenith) << ' ' << e.injectedMuonAzimuth << ' ' << e.injectedMuonEnergy << ' ';
        flux_file << cos(e.zenith) << ' ' << e.azimuth << ' ' << e.energy << ' ';
        flux_file << e.cachedConvPionWeight + e.cachedConvKaonWeight << std::endl;
    }
  }
  if(!quiet)
    std::cout << "End Loading nuSQuIDS objects." << std::endl;
  t.stop();
	if(!quiet)
		std::cout << t << std::endl;
/////////////////////////////////////////////////// VERY IMPORTANT // MUY IMPORTANTE ///////////////////////////////////////////////////
// THIS IS THE WEIGHT MAKER MASTER OBJECT // OBJECTO CREADOR DE PESADORES
	DiffuseFitWeighterMaker DFWM;
/////////////////////////////////////////////////// VERY IMPORTANT // MUY IMPORTANTE ///////////////////////////////////////////////////
  if(dump_mc_data){
    std::ofstream dump_mc = std::ofstream(output_path + "/" +"NuFSGenMC.dat");
    std::cout << "Begin Dump MC info." << std::endl;
    std::cout << "MC used " << simulation_to_load << " with DOM eff = " << std::to_string(0.9*(1.0+dataChallenge_nuisance_parameters[4])) << std::endl;
    auto weighter=DFWM(dataChallenge_nuisance_parameters);
    dump_mc << std::scientific;
    dump_mc.precision(6);
    dump_mc.fill('0');
    for(const Event& e : mainSimulation){
      if( e.energy > 4.e2 and e.energy < 2.e4 ){
        LW::Event lw_e {e.leptonEnergyFraction,
                  e.injectedEnergy,
                  e.totalColumnDepth,
                  e.inelasticityProbability,
                  e.intX,
                  e.intY,
                  e.injectedMuonEnergy,
                  e.injectedMuonZenith,
                  e.injectedMuonAzimuth,
                  static_cast<particleType>(e.primaryType),
                  e.year};
      //{
        auto w=weighter(e);
        dump_mc << static_cast<int>(e.primaryType) << ' ';
        dump_mc << e.energy << ' ' << cos(e.zenith) << ' ';
        //dump_mc << e.injectedEnergy << ' ' << e.intX << ' ' << e.intY << ' ';
        dump_mc << e.injectedEnergy << ' ' << cos(e.injectedMuonZenith) << ' ';
        //dump_mc << e.injectedMuonEnergy << ' ' << e.injectedMuonZenith << ' ';
        //dump_mc << e.inelasticityProbability << ' ' << e.totalColumnDepth << ' ';
        dump_mc << w/2. << ' ';
        //dump_mc << e.cachedConvKaonWeight << ' ';
        dump_mc << flux_pion->EvaluateFlux(lw_e) << ' ';
        dump_mc << flux_kaon->EvaluateFlux(lw_e) << std::endl;
        //dump_mc << w << std::endl;
      }
    }
    dump_mc.close();
    std::cout << "End Dump MC info." << std::endl;
    exit(0);
  }

/*
 * NOW CHRISTOPHER N WEAVER IS GOING TO CREATE HISTOGRAMS USING HIS
 * AWESOME HISTOGRAM CLASS. HERE WE CAN CHANGE THE SEARCH BINNING
 * BUT ALL IS GOOD.
 */

	using namespace phys_tools::histograms;

	std::mt19937 rng;
	if(rngSeed==0)
		rngSeed=time(0);
  rng.seed(rngSeed);
  if(!quiet)
    std::cout << "RNG seed = " << rngSeed << std::endl;

	using namespace likelihood;

	//Axes are: 
	//0: Energy Proxy (logarithmic)
	//1: cos(Zenith) (linear)
	//2: Year (linear)
	using HistType = histogram<3,entryStoringBin<std::reference_wrapper<const Event>>>;
	using phys_tools::histograms::amount;
	
	auto binner = [](HistType& h, const Event& e){
		h.add(e.energy,cos(e.zenith),e.year,amount(std::cref(e)));
	};
	
	//HistType dataHist(LogarithmicAxis(0,.1),LinearAxis(0,.1),LinearAxis(2010,1));
	//HistType dataHist(LogarithmicAxis(0,0.169897),LinearAxis(0,.05),LinearAxis(2010,1));
	HistType dataHist(LogarithmicAxis(0,0.169897),LinearAxis(0,.06),LinearAxis(2010,1));
	//HistType dataHist(LogarithmicAxis(0,0.169897),LinearAxis(0,.05),LinearAxis(2010,1));
	dataHist.getAxis(0)->setLowerLimit(minFitEnergy);
	dataHist.getAxis(0)->setUpperLimit(maxFitEnergy);
  dataHist.getAxis(1)->setLowerLimit(minCosth);
  dataHist.getAxis(1)->setUpperLimit(maxCosth);
	
  if(!quiet)
    std::cout << "Fill data histogram" << std::endl;
	bin(burnsample,dataHist,binner);

  unsigned int num_count = 0;
  for(auto event : burnsample){
    if(event.energy > 4.0e2 and event.energy < 2.0e4)
      num_count ++;
  }
  if(!quiet)
    std::cout << "after cut num cut: " << num_count << std::endl;

  /*
  // esto mejor para la casa.
  std::cout << "print data histogram" << std::endl;
  auto dho = std::ofstream("datahist");
  //for(auto bin : dataHist){
  for(auto it=dataHist.begin(); it != dataHist.end(); it++){
    auto itc=static_cast<entryStoringBin<std::reference_wrapper<const Event>>>(*it);
    dho << it.getBinEdge(0) << " " << it.getBinEdge(1)  << " " << itc.size() << std::endl;
  }
  //dho << dataHist << std::endl;
  dho.close();
  */

	bool splitSimulation=doComputeSensitivity || sampleTestStatistics;
	decltype(mainSimulation) simSubset1, simSubset2;
	if(splitSimulation){
		unsigned long idx=0;
		for(auto& item : mainSimulation){
			item.cachedNugenWeight*=2;
			item.cachedConvPionWeight*=2;
			item.cachedConvKaonWeight*=2;
			item.cachedPromptWeight*=2;
			item.cachedAstroWeight*=2;
			if((++idx)%2 == 0)
				simSubset1.push_back(item);
			else
				simSubset2.push_back(item);
		}
		mainSimulation.swap(simSubset1);
		simSubset1.clear(); //reclaim memory from duplicate data
		if(!quiet)
			std::cout << "Split simulation into 2 subsets of size " << mainSimulation.size() << " and " << simSubset2.size() << std::endl;
	}

	std::vector<HistType> allSimHists;
  {
    std::vector<std::deque<Event>*> allSimSets{&mainSimulation};
    for(const auto& simSet : allSimSets){
      allSimHists.push_back(makeEmptyHistogramCopy(dataHist));
      auto& hist=allSimHists.back();
      bin(*simSet,hist,binner);
    }
  }
  if(!quiet)
    std::cout << "End Fill simulation histogram" << std::endl;

/////////////////////////////////////////////////// VERY IMPORTANT // MUY IMPORTANTE ///////////////////////////////////////////////////
// HERE WE SPECIFY THE CONTINIOUS PARAMETER PRIORS FOR THE LIKELIHOOD // AQUI ESPECIFICAMOS LAS DISTRIBUCIONES A PRIORI DE LOS PARAMETROS CONTINUOS PARA LA FUNCION DE VEROSIMILITUD
  if(!quiet)
    std::cout << "Begin constructing priors" << std::endl;

  std::map<std::string,double> delta_alpha {
                                            {"Honda+Gaisser",8./7.},
                                            {"CombinedGHandHG_H3a_QGSJET",4./7.},
                                            {"CombinedGHandHG_H3a_SIBYLL2",8./7.},
                                            {"PolyGonato_QGSJET-II-04",0.5},
                                            {"PolyGonato_SIBYLL2",1.0},
                                            {"ZatsepinSokolskaya_pamela_QGSJET",5./7.},
                                            {"ZatsepinSokolskaya_pamela_SIBYLL2",5./7.},
                                           };

  if( delta_alpha.find(delta_model_name) == delta_alpha.end() )
    throw std::runtime_error("Jordi delta key not found. Die! Die!");

  double alpha = delta_alpha[delta_model_name];

	UniformPrior positivePrior(0.0,std::numeric_limits<double>::infinity());
	GaussianPrior normalizationPrior(1.,0.4);
	UniformPrior noPrior;
	GaussianPrior crSlopePrior(0.0,0.05);
	UniformPrior simple_domEffPrior(-.1,.3);
	GaussianPrior kaonPrior(1.0,0.1);
	GaussianPrior ZCPrior(0.0,0.038*alpha);
	GaussianPrior nanPrior(1.0,0.1);

	//auto priors=makePriorSet(positivePrior,positivePrior,positivePrior,
	//						 crSlopePrior,simple_domEffPrior,kaonPrior,positivePrior,ZCPrior);
	//auto priors=makePriorSet(positivePrior,positivePrior,positivePrior,
	//						 crSlopePrior,simple_domEffPrior,kaonPrior,nanPrior,ZCPrior);
	auto priors=makePriorSet(normalizationPrior,positivePrior,positivePrior,
							 crSlopePrior,simple_domEffPrior,kaonPrior,nanPrior,ZCPrior);
  if(!quiet)
    std::cout << "End constructing priors" << std::endl;
/////////////////////////////////////////////////// VERY IMPORTANT // MUY IMPORTANTE ///////////////////////////////////////////////////

	t.setName("");
	t.start();

	if(doDataChallenge){
    t.setName("doDataChallenge");
    t.start();
    HistType sampleHist=makeEmptyHistogramCopy(dataHist);
    std::deque<Event> DCSimulation;
    bool is_dc_simulation_already_loaded = false;

    // now lets remake the main simulation -- this already happened
    //initializeSimulationWeights(mainSimulation,PionFluxWeighter,KaonFluxWeighter,PromptFluxWeighter,osw_main);
    // make data histogram
    HistType simHist=makeEmptyHistogramCopy(dataHist);
    bin(mainSimulation,simHist,binner);

    for(double dc_seed: dataChallenge_seeds){
      sampleHist = makeEmptyHistogramCopy(dataHist);
      if(use_datachallenge_histogram){
        std::vector<Event> sample;
        //if (data_challenge_histogram_filepath == "")
        data_challenge_histogram_filepath = "dc_sample_"+std::to_string(dc_seed) + ".dat";
        std::ifstream file(data_challenge_histogram_filepath);
        if(!file.good())
          throw std::runtime_error(data_challenge_histogram_filepath + " not found.");

        double zenith,energy,year;
        while(file >> zenith >> energy >> year){
          Event e;
          e.zenith = zenith;
          e.energy = energy;
          e.year = year;
          sample.push_back(e);
        }

        bin(sample,sampleHist,binner);
        if(!quiet)
          std::cout << "Loaded DataChallenge Sample : " + data_challenge_histogram_filepath << std::endl;
      } else {
        // load auxiliary simulation from which DC will be drawn
        const bool loadTargeted=true;
        const std::vector<std::string> DCsimSetsToLoad = { dcsimulation_to_load };
        const bool is_same_simulation = (dcsimulation_to_load == simulation_to_load);
        //const bool is_same_simulation = true;
        if (is_same_simulation){
          //DCSimulation = mainSimulation;
        } else if (is_dc_simulation_already_loaded) {
          // dont do anything
        } else {
          try{
            if(readCompact){
              unsplatData(compact_data_path+"/"+dcsimulation_to_load+"_compact_data.dat",getFileChecksum(argv[0]),burnsample,DCSimulation);
            } else {
              loadSimulatedData(DCSimulation,mc_dataPath,livetime,simInfo,DCsimSetsToLoad,loadTargeted);
            }
#ifdef _FULL_EVENT_
            for(auto it = DCSimulation.begin(); it!=DCSimulation.end(); it++){
              if ((*it).azimuth < minAzimuth or (*it).azimuth > maxAzimuth)
                DCSimulation.erase(it);
              else if ((*it).cogz < minCOGZ or (*it).cogz > maxCOGZ )
                DCSimulation.erase(it);
            }
#endif
            is_dc_simulation_already_loaded = true;
          } catch(std::exception& ex){
            std::cerr << "Problem loading simulated data: " << ex.what() << std::endl;
            return(1);
          }
        }
        if(!quiet and !is_same_simulation)
          std::cout << "Loaded " << DCSimulation.size() << " events in data challenge simulation set" << std::endl;

        // put the seed into the random number generator
        rng.seed(dc_seed);
        // create null hypohtesis weighters
        std::shared_ptr<LW::SQUIDSFlux> flux_kaon_null,flux_pion_null;
        if (model_name == ""){
          flux_kaon_null = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "kaon_atmospheric_"+std::to_string(dm41sq_null)+"_"+std::to_string(th24_null)+".hdf5");
          flux_pion_null = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "pion_atmospheric_"+std::to_string(dm41sq_null)+"_"+std::to_string(th24_null)+".hdf5");
        } else {
          flux_kaon_null = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "kaon_atmospheric_"+std::to_string(dm41sq_null)+"_"+std::to_string(th24_null)+"_"+model_name+".hdf5");
          flux_pion_null = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "pion_atmospheric_"+std::to_string(dm41sq_null)+"_"+std::to_string(th24_null)+"_"+model_name+".hdf5");
        }
        //std::shared_ptr<LW::SQUIDSFlux> flux_kaon_null = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "kaon_atmospheric_"+std::to_string(dm41sq_null)+"_"+std::to_string(th24_null)+".hdf5");
        //std::shared_ptr<LW::SQUIDSFlux> flux_pion_null = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "pion_atmospheric_"+std::to_string(dm41sq_null)+"_"+std::to_string(th24_null)+".hdf5");
        // created MC weighter for DC simulation
        OversizeWeighter osw_dc(  oversize_function_dc  ) ;
        LW::mcgenWeighter dcmcw;
        for( std::string sim_name : DCsimSetsToLoad )
          dcmcw.addGenerationSpectrum(simInfo.find(sim_name)->second.details);
        // create lepton weighters for the DC simulations
        LW::LeptonWeighter PionFluxWeighter_null(flux_pion_null,xsw,dcmcw);
        LW::LeptonWeighter KaonFluxWeighter_null(flux_kaon_null,xsw,dcmcw);
        // initialize simulation weights // inicializa los pesos de la simulacion
        if (is_same_simulation)
          initializeSimulationWeights(mainSimulation,PionFluxWeighter_null,KaonFluxWeighter_null,PromptFluxWeighter,osw_dc);
        else
          initializeSimulationWeights(DCSimulation,PionFluxWeighter_null,KaonFluxWeighter_null,PromptFluxWeighter,osw_dc);
        auto weighter=DFWM(dataChallenge_nuisance_parameters);

        std::vector<double> weights;
        if (is_same_simulation)
          weights.reserve(mainSimulation.size());
        else
          weights.reserve(DCSimulation.size());
        double expected=0, expectedNuMu=0, expectedNuTau=0;
        if(is_same_simulation){
          for(const Event& e : mainSimulation){
            auto w=weighter(e);
            if(std::isnan(w) || std::isinf(w) || w<0){
              std::cout << "Bad weight!" << std::endl;
              std::cout << e.cachedConvPionWeight  << ' ' << e.cachedConvKaonWeight << ' ' << e.cachedLivetime << ' ';
              std::cout << e.energy << ' ' << e.year << ' ' << w << std::endl;
            }
            weights.push_back(w);
            //std::cout << (int) e.primaryType << std::endl;
            if(sameGeneration(particleType::NuMu,e.primaryType) or 
               (int) e.primaryType == 13 or (int) e.primaryType == -13)
              expectedNuMu+=w;
            else if(sameGeneration(particleType::NuTau,e.primaryType))
              expectedNuTau+=w;
          }
          expected= expectedNuMu + expectedNuTau;
        } else {
          for(const Event& e : DCSimulation){
            auto w=weighter(e);
            if(std::isnan(w) || std::isinf(w) || w<0){
              std::cout << "Bad weight!" << std::endl;
              std::cout << e.cachedConvPionWeight  << ' ' << e.cachedConvKaonWeight << ' ' << e.cachedLivetime << ' ';
              std::cout << e.energy << ' ' << e.year << ' ' << w << std::endl;
            }
            weights.push_back(w);
            //std::cout << (int) e.primaryType << std::endl;
            if(sameGeneration(particleType::NuMu,e.primaryType) or 
               (int) e.primaryType == 13 or (int) e.primaryType == -13)
              expectedNuMu+=w;
            else if(sameGeneration(particleType::NuTau,e.primaryType))
              expectedNuTau+=w;
          }
          expected= expectedNuMu + expectedNuTau;
        }

        if (!quiet){
          std::cout << "Expect " << expected << " events" << std::endl;
          std::cout << " of which " << expectedNuMu << " nu_mu and " << expectedNuTau << " nu_tau " << std::endl;
          if(is_same_simulation)
            std::cout << " (from weighting " << mainSimulation.size() << " simulated events)" << std::endl;
          else
            std::cout << " (from weighting " << DCSimulation.size() << " simulated events)" << std::endl;
        }

        for(auto& weight : weights)
          weight/=expected;

        std::vector<Event> sample;
        if(is_same_simulation)
          sample=likelihood::generateSample(weights,mainSimulation,expected,rng);
        else
          sample=likelihood::generateSample(weights,DCSimulation,expected,rng);

        if(!quiet)
          std::cout << " sampled " << sample.size() << " events" << std::endl;

        bin(sample,sampleHist,binner);

        if(save_dc_sample){
          //std::cout << dc_seed << std::endl;
          //if(data_challenge_histogram_filepath == "")
          data_challenge_histogram_filepath = "dc_sample_"+std::to_string(dc_seed) + ".dat";
          std::ofstream dcsample_file(data_challenge_histogram_filepath);
          for(auto& e : sample) dcsample_file << e.zenith << ' ' << e.energy << ' ' << e.year <<'\n';
          dcsample_file.close();
        }

        if(dump_data){
          const std::string sterile_params_str_null = "dm41sq_null_"+ std::to_string(dm41sq_null)+"_th24_null_"+std::to_string(th24_null);

          // i don't understand CW plotting code fully. Lets just dump all the info.
          std::ofstream dump_file(output_path + dcsimulation_to_load + "_datachallenge_data_output_"+
                                  sterile_params_str_null+"_"+model_name+"_.dat");

          for(Event& e : sample) {
            dump_file << static_cast<int>(e.primaryType) << ' ';
            dump_file << cos(e.injectedMuonZenith) << ' ' << e.injectedMuonAzimuth << ' ' << e.injectedMuonEnergy << ' ';
            dump_file << cos(e.zenith) << ' ' << e.azimuth << ' ' << e.energy << ' ';
#ifdef _FULL_EVENT_
            dump_file << e.cogx << ' ' << e.cogy << ' ' << e.cogz << ' ';
            dump_file << e.rlogl << ' ';
            dump_file << e.time << ' ';
            dump_file << e.n_chan << ' ' << e.n_chan_nodc << ' ';
            dump_file << e.intX << ' ' << e.intY << ' ';
#endif
            dump_file << 1. << std::endl;
          }

          dump_file.close();
        }
      }


      auto prob=makeLikelihoodProblem<std::reference_wrapper<const Event>,3,8>(sampleHist, {simHist}, priors, {1.0}, simpleDataWeighter(), DFWM, poissonLikelihood(), fitSeed);
      prob.setEvaluationThreadCount(evalThreads);

      std::vector<double> seed=prob.getSeed();
      std::vector<unsigned int> fixedIndices;
      for(const auto pf : fixedParams.params){
        if(!quiet)
          std::cout << "Fitting with parameter " << pf.first << " fixed to " << pf.second << std::endl;
        seed[pf.first]=pf.second;
        fixedIndices.push_back(pf.first);
      }
      fitResult fr;
      if(!use_gsl_minimizer)
        fr = doFitLBFGSB(prob,seed,fixedIndices);
      else
        fr = doFitGSL(prob,seed,fixedIndices);

      if(!quiet){
        std::cout << "Input Hypothesis: ";
        for(unsigned int i=0; i<dataChallenge_nuisance_parameters.size(); i++)
          std::cout << dataChallenge_nuisance_parameters[i] << ' ';
        std::cout << std::endl;
      }

      if(!quiet)
        std::cout << "Fitted Hypothesis: ";
      for(unsigned int i=0; i<fr.params.size(); i++)
        std::cout << fr.params[i] << ' ';
      std::cout << dm41sq << " " << th24 << " ";
      std::cout << std::setprecision(10) << fr.likelihood << std::setprecision(6) << ' ' << (fr.succeeded?"succeeded":"failed") << std::endl;

      if(dump_fit){
        const std::string sterile_params_str = "dm41sq_"+ std::to_string(dm41sq)+"_th24_"+std::to_string(th24);
        const std::string sterile_params_str_null = "dm41sq_null_"+ std::to_string(dm41sq_null)+"_th24_null_"+std::to_string(th24_null);

        // i don't understand CW plotting code fully. Lets just dump all the info.
        std::ofstream dump_file(output_path + simulation_to_load + "_simulation_output_"+
                                sterile_params_str+"_"+sterile_params_str_null+"_"+model_name+"_.dat");

        for ( Event& e : mainSimulation ) {
            auto convWeight=DFWM(fr.params);
            double weight=convWeight(e);
            dump_file << static_cast<int>(e.primaryType) << ' ';
            dump_file << cos(e.injectedMuonZenith) << ' ' << e.injectedMuonAzimuth << ' ' << e.injectedMuonEnergy << ' ';
            dump_file << cos(e.zenith) << ' ' << e.azimuth << ' ' << e.energy << ' ';
#ifdef _FULL_EVENT_
            dump_file << e.cogx << ' ' << e.cogy << ' ' << e.cogz << ' ';
            dump_file << e.rlogl << ' ';
            dump_file << e.time << ' ';
            dump_file << e.n_chan << ' ' << e.n_chan_nodc << ' ';
            dump_file << e.intX << ' ' << e.intY << ' ';
#endif
            dump_file << weight << std::endl;
        }

        dump_file.close();
      }
    }
/*
    if(!use_gsl_minimizer)
      fr = doFitLBFGSB(prob,seed,fixedIndices);
    else
      fr = doFitGSL(prob,seed,fixedIndices);

    double ldiff=fr.likelihood;
    std::cout << "Null Hypothesis: ";
    for(unsigned int i=0; i<fr.params.size(); i++)
      std::cout << fr.params[i] << ' ';
    std::cout << std::setprecision(10) << fr.likelihood << std::setprecision(6) << ' ' << (fr.succeeded?"succeeded":"failed") << std::endl;
    ldiff-=fr.likelihood;
    std::cout << "Likelihood difference: " << -ldiff << std::endl;
*/

    t.stop();
    if(!quiet)
      std::cout << t << std::endl;
    return(0);
    //data challenge ends
	}

	if(doComputeSensitivity){
		double sens=computeSensitivity(mainSimulation, simSubset2, makeEmptyHistogramCopy(dataHist), binner, priors, DFWM, fitSeed, rng, sensitivityTrials);
		std::cout << "Senstivity: " << sens << std::endl;
		return(0);
	}

  if(doBootstrap){
    double expected = burnsample.size();
    std::vector<double> weights;
    weights.reserve(burnsample.size());
    std::fill(weights.begin(),weights.end(),1.);
    auto bootstrap_sample=likelihood::generateSample(weights,burnsample,expected,rng);
    HistType sampleHist=makeEmptyHistogramCopy(dataHist);
    bin(bootstrap_sample,sampleHist,binner);

    auto prob=makeLikelihoodProblem<std::reference_wrapper<const Event>,3,8>(sampleHist, allSimHists, priors, {1.0}, simpleDataWeighter(), DFWM, poissonLikelihood(), fitSeed);
    prob.setEvaluationThreadCount(evalThreads);

    std::vector<double> seed=prob.getSeed();
    std::vector<unsigned int> fixedIndices;
    for(const auto pf : fixedParams.params){
      if(!quiet)
        std::cout << "Fitting with parameter " << pf.first << " fixed to " << pf.second << std::endl;
      seed[pf.first]=pf.second;
      fixedIndices.push_back(pf.first);
    }
    fitResult fr;
    if(!use_gsl_minimizer)
      fr = doFitLBFGSB(prob,seed,fixedIndices);
    else
      fr = doFitGSL(prob,seed,fixedIndices);

    for(unsigned int i=0; i<fr.params.size(); i++)
      std::cout << fr.params[i] << ' ';
    std::cout << dm41sq << " " << th24 << " ";
    std::cout << std::setprecision(10) << fr.likelihood << std::setprecision(6) << ' ' << (fr.succeeded?"succeeded":"failed") << std::endl;

    return 0;
  }

  t.start("FittingTime");

  if(!quiet)
    std::cout << "Making Likelihood problem. Begin fitting" << std::endl;
	auto prob=makeLikelihoodProblem<std::reference_wrapper<const Event>,3,8>(dataHist, allSimHists, priors, {1.0}, simpleDataWeighter(), DFWM, poissonLikelihood(), fitSeed);
	prob.setEvaluationThreadCount(evalThreads);

 	std::vector<double> seed=prob.getSeed();
	std::vector<unsigned int> fixedIndices;
	for(const auto pf : fixedParams.params){
		if(!quiet)
			std::cout << "Fitting with parameter " << pf.first << " fixed to " << pf.second << std::endl;
		seed[pf.first]=pf.second;
		fixedIndices.push_back(pf.first);
	}
	fitResult fr;
	if(existingBest.empty()){
    if(!use_gsl_minimizer)
      fr = doFitLBFGSB(prob,seed,fixedIndices);
    else
      fr = doFitGSL(prob,seed,fixedIndices);

		for(unsigned int i=0; i<fr.params.size(); i++)
			std::cout << fr.params[i] << ' ';
    std::cout << dm41sq << " " << th24 << " ";
		std::cout << std::setprecision(10) << fr.likelihood << std::setprecision(6) << ' ' << (fr.succeeded?"succeeded":"failed") << std::endl;

		if(!quiet){
			std::cout << fr.nEval << " likelihood evaluations" << std::endl;
			std::cout << fr.nGrad << " gradient evaluations" << std::endl;
			std::vector<double> grad=prob.gradient(fr.params);
			std::cout << "gradient = ";
			for(unsigned int i=0; i<grad.size(); i++)
				std::cout << grad[i] << ' ';
			std::cout << std::endl;
		}
	}
	else{
		if(!quiet)
			std::cout << "Using existing fit: " << existingBest_ << std::endl;
		fr.params=existingBest;
		fr.likelihood=-prob.evaluateLikelihood(fr.params);
		std::cout << "Likelihood: " << fr.likelihood << std::endl;
	}

  t.stop();
  if(!quiet)
    std::cout << t << std::endl;

  // Paramos capitan. Paramos.
  //exit(0);

  {
    std::cout << "print mc histogram" << std::endl;
    auto mho = std::ofstream("mchist");
    //for(auto bin : dataHist){
    auto convWeight=DFWM(fr.params);
    for(auto it=allSimHists[0].begin(); it != allSimHists[0].end(); it++){
      auto itc=static_cast<entryStoringBin<std::reference_wrapper<const Event>>>(*it);
      mho << it.getBinEdge(0) << " " << it.getBinEdge(1)  << " ";
      double sum = 0.;
      for(auto event : itc.entries()){
          sum += convWeight(event);
      }
      mho << sum << std::endl;
    }
    mho.close();
  }

/*
	//check GoF:
	if(!quiet){
		double ref=saturatedPoisson(prob,fr.params);
		std::cout << "Reference: " << ref << std::endl;
		//unsigned int dof=5;
		unsigned int dof=prob.likelihoodTerms()-(7-fixedIndices.size());
		std::cout << " Problem has " << dof << " degrees of freedom" << std::endl;
		double R=2*(fr.likelihood-ref);
		std::cout << " R=" << R << std::endl;
		auto dist=boost::math::chi_squared_distribution<double>(dof*1.3);
		double intProb=boost::math::cdf(dist,R);
		std::cout << " P=" << 1-intProb << std::endl;
		
		std::cout << "----\n";
		auto sChi2prob=prob.makeAlternateLikelihood(likelihood::chi2Likelihood());
		double chi2=sChi2prob.evaluateLikelihood(fr.params);
		std::cout << "Chi^2: " << chi2 << std::endl;
		auto distTheo=boost::math::chi_squared_distribution<double>(dof);
		intProb=boost::math::cdf(dist,chi2);
		std::cout << " P=" << 1-intProb << std::endl;
		
		auto dprob=prob.makeAlternateLikelihood(likelihood::dimaLikelihood());
		double dTS=dprob.evaluateLikelihood(fr.params,false);
		std::cout << "  DimaTS: " << dTS << std::endl;
	}
*/

	if(computeWilksSig)
		wilksNonZeroSig(prob,fr.params,fixedIndices);
	
	if(sampleTestStatistics){
		std::vector<double> nullHyp=fr.params;
		nullHyp[2]=0.0;
		FCNonZeroSig(mainSimulation, simSubset2, dataHist, binner, priors, DFWM, fitSeed, rng, nullHyp, burnsample.size(), sampleTestStatistics, fixedParams);
	}

	if(computeParameterErrors){ //find error intervals the fast way
    t.setName("computeParameterErrors");
    t.start();
		auto makeScanFunc=[&](unsigned int idx){
			auto seed=fr.params;
			return([&,seed,idx](double paramValue) mutable{
				seed[idx]=paramValue;
				auto fi=fixedIndices;
				fi.push_back(idx);
        fitResult fr;
        if(!use_gsl_minimizer)
          fr = doFitLBFGSB(prob,seed,fi);
        else
          fr = doFitGSL(prob,seed,fi);
				return(fr.likelihood);
			});
		};
		
		auto findParameterError=[&](unsigned int idx, double guessLow, double guessHigh, double precision, double limitLow, double limitHigh){
			auto astroScanFunc=makeScanFunc(idx);
			double middle=fr.params[idx];
			
			std::pair<double,double> domain=bracketRoot(astroScanFunc, fr.likelihood+0.5, middle, guessLow, fr.likelihood, limitLow);
			double lower;
			if(domain.second<fr.likelihood+0.5)
				lower=domain.first;
			else
				lower=findRoot(astroScanFunc, fr.likelihood+0.5, domain.first, middle, domain.second, fr.likelihood, precision, true);
			
			domain=bracketRoot(astroScanFunc, fr.likelihood+0.5, middle, guessHigh, fr.likelihood, limitHigh);
			double upper;
			if(domain.second<fr.likelihood+0.5)
				upper=domain.first;
			else
				upper=findRoot(astroScanFunc, fr.likelihood+0.5, middle, domain.first, fr.likelihood, domain.second, precision, true);
			
			std::cout << "Error interval for parameter " << idx << " is [" << lower << ',' << upper << ']' << std::endl;
		};
		
		if(std::find(fixedIndices.begin(),fixedIndices.end(),0)==fixedIndices.end())
			findParameterError(0,fr.params[0]-.2,fr.params[0]+.2,5e-3,0,std::numeric_limits<double>::infinity());
    /*
		if(std::find(fixedIndices.begin(),fixedIndices.end(),1)==fixedIndices.end())
			findParameterError(1,0,fr.params[1]*2,1e-2,0,std::numeric_limits<double>::infinity());
		if(std::find(fixedIndices.begin(),fixedIndices.end(),2)==fixedIndices.end())
			findParameterError(2,fr.params[2]/2,fr.params[2]*2,1e-2,0,std::numeric_limits<double>::infinity());
    */
		if(std::find(fixedIndices.begin(),fixedIndices.end(),3)==fixedIndices.end())
			findParameterError(3,fr.params[3]-.01,fr.params[3]+.01,1e-3,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());
		if(std::find(fixedIndices.begin(),fixedIndices.end(),4)==fixedIndices.end())
			findParameterError(4,fr.params[4]-.005,fr.params[4]+.005,5e-4,-.1,.3);
		if(std::find(fixedIndices.begin(),fixedIndices.end(),5)==fixedIndices.end())
			findParameterError(5,fr.params[5]-.1,fr.params[5]+.1,1e-2,0,std::numeric_limits<double>::infinity());
		if(std::find(fixedIndices.begin(),fixedIndices.end(),6)==fixedIndices.end())
			findParameterError(6,fr.params[6]-.1,fr.params[6]+.1,1e-2,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());
		if(std::find(fixedIndices.begin(),fixedIndices.end(),7)==fixedIndices.end())
			findParameterError(7,fr.params[7]-.1,fr.params[7]+.1,1e-2,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());

    t.stop();
    if(!quiet)
      std::cout << t << std::endl;
	}

/*
  std::cout << "Scanning DOM efficiency" << std::endl;
  seed=fr.params;
  for(double didx=-0.2; didx<=.2; didx+=.01){
      fr.params={1,1,1,0,didx,1,1};
      fr.likelihood=-prob.evaluateLikelihood(fr.params);
      //std::cout << "Likelihood: " << fr.likelihood << std::endl;
      //seed[6]=didx;
      //fitResult fr = doFitLBFGSB(prob,seed,{6});
      for(unsigned int  i=0; i<fr.params.size(); i++)
          std::cout << fr.params[i] << ' ';
      std::cout << std::setprecision(10) << fr.likelihood << ' ' << std::setprecision(6) << ' ' << (fr.succeeded?"succeeded":"failed") << " ";//std::endl;
      std::vector<double> grad = prob.gradient(fr.params);
			std::cout << "gradient = ";
			for(unsigned int i=0; i<grad.size(); i++)
				std::cout << grad[i] << ' ';
			std::cout << std::endl;
      //seed=fr.params;
  }
*/
/*
  std::cout << "Scanning neutrino/antineutrino ratio" << std::endl;
  seed=fr.params;
  for(double didx=.5; didx<=1.5; didx+=.025){
      fr.params={1,1,1,0,.19,1,didx};
      fr.likelihood=-prob.evaluateLikelihood(fr.params);
      //std::cout << "Likelihood: " << fr.likelihood << std::endl;
      //seed[6]=didx;
      //fitResult fr = doFitLBFGSB(prob,seed,{6});
      for(unsigned int  i=0; i<fr.params.size(); i++)
          std::cout << fr.params[i] << ' ';
      std::cout << std::setprecision(10) << fr.likelihood << ' ' << std::setprecision(6) << ' ' << (fr.succeeded?"succeeded":"failed") << " ";//std::endl;
      std::vector<double> grad = prob.gradient(fr.params);
			std::cout << "gradient = ";
			for(unsigned int i=0; i<grad.size(); i++)
				std::cout << grad[i] << ' ';
			std::cout << std::endl;
      //seed=fr.params;
  }
*/

/*====================================================== CARLOS DUMP  ================================================================*/
/*====================================================== CARLOS DUMP  ================================================================*/
/*====================================================== CARLOS DUMP  ================================================================*/

t.start("dumpSimulationOutput");

const std::string sterile_params_str = "dm41sq_"+ std::to_string(dm41sq)+"_th24_"+std::to_string(th24);

// i don't understand CW plotting code fully. Lets just dump all the info.
std::ofstream dump_file(output_path + simulation_to_load + "_simulation_output_"+
                        sterile_params_str+"_"+model_name+"_.dat");

for ( Event& e : mainSimulation ) {
		auto convWeight=DFWM(fr.params);
    double weight=convWeight(e);
    //dump_file << static_cast<int>(e.primaryType) << ' ';
    //dump_file << cos(e.injectedMuonZenith) << ' ' << e.injectedMuonAzimuth << ' ' << e.injectedMuonEnergy << ' ';
    //dump_file << cos(e.zenith) << ' ' << e.azimuth << ' ' << e.energy << ' ';
    dump_file << cos(e.zenith) << ' ' << e.energy << ' ';
#ifdef _FULL_EVENT_
    dump_file << e.cogx << ' ' << e.cogy << ' ' << e.cogz << ' ';
    dump_file << e.rlogl << ' ';
    dump_file << e.time << ' ';
    dump_file << e.n_chan << ' ' << e.n_chan_nodc << ' ';
    dump_file << e.intX << ' ' << e.intY << ' ';
#endif
    dump_file << weight << std::endl;
}

dump_file.close();

// save best fit parameters to file
std::ofstream results_file(output_path + simulation_to_load + "_fit_result_"+
                           sterile_params_str+"_"+model_name+"_.dat");

for(unsigned int i=0; i<fr.params.size(); i++)
  results_file << fr.params[i] << ' ';
results_file << std::setprecision(10) << fr.likelihood << std::setprecision(6) << std::endl;

results_file.close();

t.stop();
if(!quiet)
  std::cout << t << std::endl;

/*====================================================== CWEAVER PLOTTING SCRIPTS ================================================================*/
/*====================================================== CWEAVER PLOTTING SCRIPTS ================================================================*/
/*====================================================== CWEAVER PLOTTING SCRIPTS ================================================================*/

	if(drawPlots){
		
    // plotting cos ranges
		double cosZenMin=minCosth, cosZenMax=maxCosth;
    bool plot_prompt = true;
	
		using namespace phys_tools;
		using namespace phys_tools::gnuplot;
		LinearAxis zenAxis(0.0,.02); // zenith energy binning for histogram
		LinearAxis aziAxis(0.0,pi<double>()/10.); // azimuth binning for histogram
		LogarithmicAxis enAxis(3,.1);
		enAxis.setLowerLimit(1e2);
		enAxis.setUpperLimit(1e7);
		
		LinearAxis coarseZenAxis(0.0,.1);
		
		LogarithmicAxis mcEnAxis(3,.1);
		mcEnAxis.setLowerLimit(1e1);
		mcEnAxis.setUpperLimit(1e8);
		
    FilledErrorHist<sqErrorValue> convAzimuth("Conventional atmospheric",Format::FILLED_CURVES,aziAxis);

		Plottable1DHistogram sumZenith("Sum of predictions",Format::LEFT_STEPS,zenAxis);
		FilledErrorHist<sqErrorValue> convZenith("Conventional atmospheric",Format::FILLED_CURVES,zenAxis),
		                              promptZenith("Prompt atmospheric",Format::FILLED_CURVES,zenAxis),
		                              astroZenith("E^{-2} Astrophysical",Format::FILLED_CURVES,zenAxis);
		FilledErrorHist<fcErrorValue> dataZenith("Experimental data",Format::POINTS,zenAxis);
		Plottable1DHistogram sumEnergy("Sum of predictions",Format::LEFT_STEPS,enAxis);
		FilledErrorHist<sqErrorValue> convEnergy("Conventional atmospheric",Format::FILLED_CURVES,enAxis),
		                              promptEnergy("Prompt atmospheric",Format::FILLED_CURVES,enAxis),
		                              astroEnergy("E^{-2} Astrophysical",Format::FILLED_CURVES,enAxis);
		FilledErrorHist<fcErrorValue> dataEnergy("Experimental data",Format::POINTS,enAxis);
		FilledErrorHist<fcErrorValue> dataAzimuth("Experimental data",Format::POINTS,aziAxis);
		                     
		Plottable1DHistogram convTrueEnergy("Honda+GaisserH3a Knee",Format::LEFT_STEPS,mcEnAxis),
		                     promptTrueEnergy("Enberg+GaisserH3a Knee",Format::LEFT_STEPS,mcEnAxis);
		//Plottable1DHistogram convMuonEnergy("Honda+GaisserH3a Knee",Format::LEFT_STEPS,mcEnAxis), 
		//                     promptMuonEnergy("Enberg+GaisserH3a Knee",Format::LEFT_STEPS,mcEnAxis);
		Plottable2DHistogram coarseDataZenEn("",enAxis,coarseZenAxis);
		Plottable2DHistogram coarseSimZenEn("",enAxis,coarseZenAxis);
		
		Plottable2DHistogram dataZenEn("",enAxis,zenAxis);
		Plottable2DHistogram convZenEn("",enAxis,zenAxis);
		Plottable2DHistogram promptZenEn("",enAxis,zenAxis);
		Plottable2DHistogram astroZenEn("",enAxis,zenAxis);
		dataZenEn.setUseContentScaling(false);
		convZenEn.setUseContentScaling(false);
		promptZenEn.setUseContentScaling(false);
		astroZenEn.setUseContentScaling(false);
		
		//dataZenith.format().error_style(Format::ERROR_BARS);
		//dataEnergy.format().error_style(Format::ERROR_BARS);
		
		//minFitEnergy=3.5e4;
		
		double tmpDataRate=0.0, tmpSimRate=0.0;
		if(!quiet)
			std::cout << "binning data" << std::endl;
		for(const Event& e : burnsample){
		//for(const Event& e : sample){
			if(e.energy>=minFitEnergy && e.energy<=maxFitEnergy)
				dataZenith.add(cos(e.zenith));
			if(cos(e.zenith)>=cosZenMin && cos(e.zenith)<=cosZenMax)
				dataEnergy.add(e.energy);
      dataAzimuth.add(e.azimuth);
			coarseDataZenEn.add(e.energy,cos(e.zenith));
			dataZenEn.add(e.energy,cos(e.zenith));
			if(e.energy>7e2 && e.energy<8e2)
				tmpDataRate+=1;
		}
		/*#warning treating simulation sample as true data for plotting
		for(const Event& e : mainSimulation){
			if(e.energy>=minFitEnergy && e.energy<=maxFitEnergy)
				dataZenith.add(cos(e.zenith),PlottableHistogramFCErr::amount(weighter(e)));
			if(cos(e.zenith)>=cosZenMin && cos(e.zenith)<=cosZenMax)
				dataEnergy.add(e.energy,PlottableHistogramFCErr::amount(weighter(e)));
			dataZenEn.add(e.energy,cos(e.zenith),Plottable2DHistogram::amount(weighter(e)));
		}*/
		
		auto convExp=fr.params;
		convExp[1]=0; convExp[2]=0;
		auto convWeight=DFWM(convExp);
		auto promptExp=fr.params;
		promptExp[0]=0; promptExp[2]=0;
		auto promptWeight=DFWM(promptExp);
		auto sigExp=fr.params;
		sigExp[0]=0; sigExp[1]=0;
		auto signalWeight=DFWM(sigExp);
		auto sumWeight=DFWM(fr.params);
		
		//histogram<2> muonVsNeutrinoEnergy(LogarithmicAxis(0,.1),LogarithmicAxis(0,.1));
		//muonVsNeutrinoEnergy.getAxis(0)->setLowerLimit(1e2);
		
		if(!quiet)
			std::cout << "binning simulation" << std::endl;
		for(const Event& e : mainSimulation){
			if(e.energy>maxFitEnergy)
				continue;
			
			if(e.energy<1e2 || e.energy>1e8)
				continue;
			
			double c=convWeight(e);
			double p=promptWeight(e);
			double a=signalWeight(e);
			
			if(std::isnan(c)){
				std::cout << "bad event: " << e.injectedMuonZenith << ' ' << e.injectedMuonEnergy << ' ' << e.zenith << ' ' << e.energy << std::endl;
        std::cout << e.cachedConvPionWeight << ' ' << e.cachedPromptWeight << ' ' << e.cachedAstroWeight << std::endl;
				std::cout << ' ' << e.cachedConvDOMEff.baseRate << ' ' << e.cachedPromptDOMEff.baseRate << ' ' << e.cachedAstroDOMEff.baseRate << std::endl;
				throw std::runtime_error("bad event weight");
			}
			if(e.energy>=minFitEnergy && e.energy<=maxFitEnergy){
				convZenith.add(cos(e.zenith),amount(c));
				promptZenith.add(cos(e.zenith),amount(p));
				astroZenith.add(cos(e.zenith),amount(a));
				sumZenith.add(cos(e.zenith),amount(c+promptWeight(e)+signalWeight(e)));
			}
			if(cos(e.zenith)>=cosZenMin && cos(e.zenith)<=cosZenMax){
				convEnergy.add(e.energy,amount(c));
				promptEnergy.add(e.energy,amount(p));
				astroEnergy.add(e.energy,amount(a));
				sumEnergy.add(e.energy,amount(c+promptWeight(e)+signalWeight(e)));
			}
        convAzimuth.add(e.azimuth,amount(c));

			
			coarseSimZenEn.add(e.energy,cos(e.zenith),amount(c+p+a));
			convZenEn.add(e.energy,cos(e.zenith),amount(c));
			promptZenEn.add(e.energy,cos(e.zenith),amount(p));
			astroZenEn.add(e.energy,cos(e.zenith),amount(a));
			
			convTrueEnergy.add(e.injectedEnergy,amount(c));
			promptTrueEnergy.add(e.injectedEnergy,amount(promptWeight(e)));
			//convMuonEnergy.add(e.muonEntryEnergy,amount(c));
			//promptMuonEnergy.add(e.muonEntryEnergy,amount(promptWeight(e)));
			
			//muonVsNeutrinoEnergy.add(e.muonEntryEnergy,e.primaryEnergy,amount(signalWeight(e)));
			if(e.energy>7e2 && e.energy<8e2)
				tmpSimRate+=c+p+a;
		}
		
		{
			auto findMinConf=[](double obs, double exp)->double{
				double min=0, max=1, test;
				while((max-min)>1e-4){
					test=(min+max)/2;
					interval<unsigned int> accept=fcPoissonAcceptanceInterval(exp, test, 0);
					(accept.contains(obs) ? max : min)=test;
				}
				return(test);
			};
			dataEnergy.setUseContentScaling(false);
			sumEnergy.setUseContentScaling(false);
			boost::math::normal_distribution<double> gauss;
			for(auto it=dataEnergy.cbegin(), end=dataEnergy.cend(); it!=end; it++){
				auto obs=*it;
				auto expIt=sumEnergy.findBinIterator(it);
				double exp=(expIt!=sumEnergy.end()?*expIt:0.0);
				double conf=findMinConf(obs,exp);
				std::cout << it.getBinEdge(0) << ' ' << obs << ' ' << exp << ' ' << conf << ' ' << boost::math::quantile(gauss,(1+conf)/2) << std::endl;
			}
		}
		
		if(!quiet)
			std::cout << "making plots" << std::endl;
		Gnuplot g;
		g.set_terminal(phys_tools::gnuplot::Terminal("postscript eps enhanced color solid lw 2 font \"Helvetica\" 20"));
		auto black=rgb_color("black");
		auto red=rgb_color("red");
		auto green=rgb_color(0x00CC00);
		auto blue=rgb_color("blue");
		auto purple=rgb_color("purple");
		auto royalblue=rgb_color("royalblue");
		auto forestgreen=rgb_color("forest-green");
		auto darkred=rgb_color(0xC00000);
		
		g.line_style(2).color(black);
		g.set_bar_size(0);
		g.x_axis().minor_tics(10);
		g.y_axis().label("Events").minor_tics(10);
		
		dataAzimuth.format().error_style(Format::Y_ERROR_BARS).line_color(black).line_width(1).point_type(7).point_size(1.5);
		dataZenith.format().error_style(Format::Y_ERROR_BARS).line_color(black).line_width(1).point_type(7).point_size(1.5);
		dataEnergy.format().error_style(Format::Y_ERROR_BARS).line_color(black).line_width(1).point_type(7).point_size(1.5);
		convEnergy.setUseContentScaling(false);

    if(!quiet)
      std::cout << "convEnergy.integral() = " << convEnergy.integral() << std::endl;

		promptEnergy.setUseContentScaling(false);
		astroEnergy.setUseContentScaling(false);
		sumEnergy.format().line_color(purple).line_width(4);
		sumEnergy.setUseContentScaling(false);
		sumZenith.format().line_color(purple).line_width(4);
		sumZenith.setUseContentScaling(false);
		
		convAzimuth.format().fill_color(red).line_color(darkred).line_width(1);
		convZenith.format().fill_color(red).line_color(darkred).line_width(1);
		promptZenith.format().fill_color(royalblue).line_color(blue).line_width(1);
		astroZenith.format().fill_color(green).line_color(forestgreen).line_width(1);
		
		dataEnergy.setUseContentScaling(false);
		
		convEnergy.format().fill_color(red).line_color(darkred).line_width(1);
		promptEnergy.format().fill_color(royalblue).line_color(blue).line_width(1);
		astroEnergy.format().fill_color(green).line_color(forestgreen).line_width(1);
		
		coarseDataZenEn.setUseContentScaling(false);
		coarseSimZenEn.setUseContentScaling(false);
		//std::cout << '[' << dataZenEn.range(0).first << ',' << dataZenEn.range(0).second << "] x [" << dataZenEn.range(1).first << ',' << dataZenEn.range(1).second << ']' << std::endl;
		//std::cout << '[' << simZenEn.range(0).first << ',' << simZenEn.range(0).second << "] x [" << simZenEn.range(1).first << ',' << simZenEn.range(1).second << ']' << std::endl;
		//std::cout << "dataZenEn: \n" << dataZenEn << std::endl;
		//std::cout << "simZenEn: \n" << simZenEn << std::endl;
		
		convTrueEnergy.format().line_color(red).line_width(2);
		promptTrueEnergy.format().line_color(royalblue).line_width(2);

    // plot azimuth distribution
		g.manual_settings("set key left Left reverse top");
		g.y_axis().label("Events");
		g.x_axis().label("Azimuth [rad]");
    g.plot(convAzimuth,dataAzimuth,plot_path+"azimuth_"+sterile_params_str+"_"+model_name+".eps");
	  if(!quiet)	
      std::cout << "zenith.eps" << std::endl;
		g.manual_settings("set key left Left reverse top");
		//g.label(1).text("IceCube Preliminary").at("first -0.9, first 200").text_color(red).font("Monaco",20);
		g.x_axis().label("cos(Reconstructed zenith angle)");
		g.y_axis().label("Events");
		//g.y_axis().label("Events").logscale(true).auto_max();
		{//For APS slide
			//dataZenith.title("'Burnsample' Experimental Data");
			//dataZenith.title("IC79+IC86 simulated data realization");
			//convZenith.title("Conventional atmospheric prediction");
			//convZenith.format().line_width(3);
		}
		//g.plot(convZenith, promptZenith, astroZenith/*, sumZenith*/, dataZenith, plot_path+"zenith.eps");
    //if(plot_prompt)
    //  g.plot(convZenith, promptZenith,sumZenith, dataZenith, plot_path+"zenith_"+sterile_params_str+"_"+model_name+".eps");
    //else
    g.plot(convZenith/*, sumZenith*/, dataZenith, plot_path+"zenith_"+sterile_params_str+"_"+model_name+".eps");
	  if(!quiet)	
      std::cout << "zenith_ratio.eps" << std::endl;
    auto zenith_ratio = convZenith.getCenter()/dataZenith.getCenter();
		g.manual_settings("set key left Left reverse top");
		g.x_axis().label("cos(Reconstructed zenith angle)");
		g.y_axis().label("data/MC");
    g.plot(Plottable1DHistogram("",Format::LEFT_STEPS,zenith_ratio),plot_path+"zenith_ratio_"+sterile_params_str+"_"+model_name+".eps");

		std::cout << "zenith_log.eps" << std::endl;
		g.manual_settings("set key left Left reverse top invert");
		g.x_axis().label("cos(Reconstructed zenith angle)");
		g.y_axis().label("Events");
		g.y_axis().logscale(true).format("10^{%L}").max(2e4);
		g.label(1).at(Coordinate(-0.9,0.4));
		//g.plot(convZenith, promptZenith, astroZenith, sumZenith, dataZenith, plot_path+"zenith_log_"+sterile_params_str+".eps");
    if(plot_prompt)
      g.plot(convZenith, promptZenith, sumZenith, dataZenith, plot_path+"zenith_log_"+sterile_params_str+"_"+model_name+".eps");
    else
      g.plot(convZenith, dataZenith, plot_path+"zenith_log_"+sterile_params_str+"_"+model_name+".eps");
		g.unset_label(1);
		
		//factor 18.532 (19.649)
		//g.manual_settings("set bars small \n unset link y \n set y2range [.19649:19649] \n set y2tics \n set format y2 '10^{%L}' \n set y2label 'Events, scaled to 2 years' \n set ytics nomirror");
		//g.x_axis().label("MuEx Energy (arb. units)").logscale(true).format("10^{%L}").min(2e2);
		//g.x_axis().label("MuEx Energy (GeV)").logscale(true).format("10^{%L}").min(2e2).max(maxFitEnergy);
		g.y_axis().label("Events").logscale(true).format("10^{%L}").max(1e4).min(.01);
		//g.y_axis().label("Events").logscale(false).format("%g").max(7000).min(0);
		g.x_axis().label("MuEx Energy (GeV)").logscale(true).format("10^{%L}").min(minFitEnergy).max(1.0e6);
		//g.y_axis().label("Events, 805 hours");
		{//For APS slide
			convEnergy.title("Conventional atmospheric");
			promptEnergy.title("Prompt atmospheric");
			astroEnergy.title("E^{-2} astrophysical");
			g.x_axis().label("Muon Energy Proxy (arb. units)");
		}
		g.plot(convEnergy,promptEnergy,astroEnergy,plot_path+"energy_no_data_"+sterile_params_str+"_"+model_name+".eps");
		
		g.manual_settings("set key right top invert");
		//g.label(1).text("IceCube Preliminary").at("first 1e3, first .1").text_color(red).font("Monaco",20);
		//g.plot(convEnergy, promptEnergy, astroEnergy, sumEnergy, dataEnergy, plot_path+"energy_"+sterile_params_str+".eps");
    if(plot_prompt)
      g.plot(convEnergy, promptEnergy, sumEnergy, dataEnergy, plot_path+"energy_"+sterile_params_str+"_"+model_name+".eps");
    else
      g.plot(convEnergy, dataEnergy, plot_path+"energy_"+sterile_params_str+"_"+model_name+".eps");
		g.unset_label(1);

    auto energy_ratio = convEnergy.getCenter()/dataEnergy.getCenter();
		g.y_axis().label("data/MC").logscale(false).format("%g").auto_min().auto_max();
    g.plot(Plottable1DHistogram("",Format::LEFT_STEPS,energy_ratio),plot_path+"energy_ratio_"+sterile_params_str+"_"+model_name+".eps");

		
		auto cumulativeConvEnergy=convEnergy;
		accumulateHistogramReverse(cumulativeConvEnergy);
		auto cumulativeAstroEnergy=astroEnergy;
		accumulateHistogramReverse(cumulativeAstroEnergy);
		auto cumulativeSumEnergy=sumEnergy;
		accumulateHistogramReverse(cumulativeSumEnergy);
		
		/*g.manual_settings("unset label 1 \n set key right top \n set bars small \n set mxtics 10 \n set mytics 10");
		Plottable1DHistogram backgroundFraction("Conventional Atmospheric Fraction",Format::LEFT_STEPS,(const histogram<1>&)cumulativeConvEnergy/(const histogram<1>&)cumulativeSumEnergy);
		Plottable1DHistogram signalFraction("Astrophysical Fraction",Format::LEFT_STEPS,(const histogram<1>&)cumulativeAstroEnergy/(const histogram<1>&)cumulativeSumEnergy);
		backgroundFraction.format().line_width(2);
		signalFraction.format().line_width(2);
		//g.y_axis().label("Fraction of Fitted Spectrum").logscale(true).format("10^{%L}").max().min(.01);
		g.y_axis().label("Fraction of Fitted Spectrum").logscale(false).format("%g").max(1.2).min(0);
		g.plot(backgroundFraction,signalFraction,"signal_ratio_energy.eps");*/
		
		g.manual_settings("set pm3d map corners2color c1 \n set palette defined ( 0 \"white\", .01 \"blue\", 3 \"green\", 6 \"yellow\", 10 \"red\")");
		{
			g.x_axis().label("Muon Energy (GeV)").logscale(true).format("10^{%L}").min(minFitEnergy);
			g.y_axis().label("cos(Zenith)").logscale(false).format("%g").auto_min().auto_max();
			g.cb_axis().logscale(true).min(1e-2).max(1e3);
			//g.cb_axis().logscale(true).min(1e-2).max(5e2);
			g.title("Experimental Distribution");
			g.plot3d(dataZenEn,plot_path+"data_en_zen_"+sterile_params_str+"_"+model_name+".eps");
			//std::cout << "dataZenEn: \n" << dataZenEn << std::endl;
			g.plot3d(coarseDataZenEn,plot_path+"data_en_zen_coarse_"+sterile_params_str+"_"+model_name+".eps");
			
			g.title("Simulated Conventional Atmospheric Distribution");
			g.plot3d(convZenEn,plot_path+"conv_en_zen.eps");
			{
				std::ofstream outfile(output_path+simulation_to_load+"_conv_en_zen_"+sterile_params_str+"_"+model_name+".txt");
				outfile << convZenEn;
			}

      /*
			g.title("Simulated Prompt Atmospheric Distribution");
			g.plot3d(promptZenEn,plot_path+"prompt_en_zen.eps");
			{
				std::ofstream outfile(output_path+"prompt_en_zen.txt");
				outfile << promptZenEn;
			}
			g.title("Simulated Astrophysical Distribution");
			g.plot3d(astroZenEn,plot_path+"astro_en_zen.eps");
			{
				std::ofstream outfile(output_path+"astro_en_zen.txt");
				outfile << astroZenEn;
			}
      */
			g.title("Simulated Distribution");
			//std::cout << "simZenEn: \n" << simZenEn << std::endl;
			g.plot3d(coarseSimZenEn,plot_path+"sim_en_zen_coarse_"+sterile_params_str+"_"+model_name+".eps");
			g.cb_axis().logscale(false).auto_min().auto_max();
			g.title("Ratio of Data to Simulation");
			//g.cb_axis().logscale(false).min(0).max(3);
			g.cb_axis().logscale(true).min(1).max(4);
			Plottable2DHistogram dataSimZenEnRatio("",(((const histogram<2>&)coarseDataZenEn)/((const histogram<2>&)coarseSimZenEn)));
			//dataSimZenEnRatio.setUseContentScaling(false); //For reasons that need to be debugged, using this would break everything
			//std::cout << "dataSimZenEnRatio: \n" << dataSimZenEnRatio << std::endl;
			g.plot3d(dataSimZenEnRatio,plot_path+"data_sim_ratio_coarse_"+sterile_params_str+"_"+model_name+".eps");
			g.title("");
			
			std::cout << "evaluating binwise likelihood" << std::endl;
			auto sol=fr.params;
			sol.push_back(0);
			std::pair<phys_tools::histograms::histogram<3,double>,double> likelihoodBits
			=prob.evaluateLikelihoodHistogram<double>(sol);
			//=prob.evaluateLikelihoodHistogram<double>(std::vector<double>{finalParams.Value(0),0,0,finalParams.Value(3),finalParams.Value(4),finalParams.Value(5)});
			likelihoodBits.first.setUseContentScaling(false);
			//std::cout << "Likelihood map: \n" << likelihoodBits.first << std::endl;
			//std::cout << " Sum: " << likelihoodBits.first.integral() << std::endl;
			//std::cout << " overflow/underflow: " << likelihoodBits.second << std::endl; //LIES
			//std::cout << "Complete likelihood: " << prob.evaluateLikelihood<double>(fr.params) << std::endl;
			
			//g.x_axis().min(likelihoodBits.first.getBinEdge(0,0)).auto_max();
			g.x_axis().min(likelihoodBits.first.getBinEdge(0,0)).max(maxFitEnergy);
			g.y_axis().label("cos(Zenith)").logscale(false).format("%g").min(-1).max(0.2);
			g.cb_axis().logscale(false).auto_min().auto_max();
			Plottable2DHistogram binwise_likelihood_map("",likelihoodBits.first.project(2));
			binwise_likelihood_map.setUseContentScaling(false);
			g.plot3d(binwise_likelihood_map,plot_path+"binwise_likelihood_map_"+sterile_params_str+"_"+model_name+".eps");
			std::cout << "Total likeihood: " << likelihoodBits.first.integral() << std::endl;
			
			std::cout << "evaluating binwise saturated likelihood" << std::endl;
			auto spProb=prob.makeAlternateLikelihood(likelihood::saturatedPoissonLikelihood());
			std::pair<phys_tools::histograms::histogram<3,double>,double> refLikelihoodBits
			=spProb.evaluateLikelihoodHistogram<double>(sol);
			
			std::cout << "Saturated likelihood: " << refLikelihoodBits.first.integral() << std::endl;
			
			Plottable2DHistogram binwise_likelihood_sat_map("",refLikelihoodBits.first.project(2)-likelihoodBits.first.project(2));
			binwise_likelihood_sat_map.setUseContentScaling(false);
			g.plot3d(binwise_likelihood_sat_map,plot_path+"binwise_likelihood_sat_map_"+sterile_params_str+"_"+model_name+".eps");
			
			histogram<3> llhSatDiff(refLikelihoodBits.first-likelihoodBits.first);
			Plottable1DHistogram diffContributions("",Format::LEFT_STEPS,LinearAxis(0,.25));
			for(const auto& diff : llhSatDiff)
				diffContributions.add(diff);
			g.manual_settings("");
			g.x_axis().label("Contribution to likelihood difference").logscale(false).format("%g").min(0).auto_max();
			g.y_axis().label("Number of bins").logscale(true).format("10^{%L}").min(0.1).auto_max();
			g.plot(diffContributions,plot_path+"binwise_likelihood_sat_map_dist_"+sterile_params_str+"_"+model_name+".eps");
		g.manual_settings("");
    }
		
		//VerticalLinePlot medianConvEn(2020.7,"Median Conv. Energy: 2.02 Tev");
		//medianConvEn.format().line_color(forestgreen).line_width(2);
		//VerticalLinePlot medianPromptEn(7887,"Median Prompt Energy: 7.887 Tev");
		//medianPromptEn.format().line_color(royalblue).line_width(2);
		//g.title("MC True Neutrino Energy");
		//g.x_axis().label("Neutrino Energy (GeV)").min(1e2);
		//g.y_axis().logscale(true).format("10^{%L}").auto_min().auto_max();
		//g.plot(convTrueEnergy,medianConvEn,promptTrueEnergy,medianPromptEn,plot_path+"mc_energy.eps");
		
		//VerticalLinePlot medianConvMuonEn(791.876,"Median Conv. Muon Energy: 792 GeV");
		//medianConvMuonEn.format().line_color("rgb \"forest-green\"").line_width(2);
		//g.title("MC True Muon Energy (at detector entry)");
		//g.x_axis().label("Muon Energy (GeV)");
		//g.plot(convMuonEnergy,medianConvMuonEn,promptMuonEnergy,"mc_muon_energy.eps");
		
		//normalizeSlices(muonVsNeutrinoEnergy);
		//g.x_axis().logscale(true).label("Muon Energy at detector entry (GeV)").format("10^{%L}").auto_min().auto_max();
		//g.y_axis().logscale(true).label("Neutrino Energy (GeV)").format("10^{%L}").auto_min().auto_max();
		//g.cb_axis().logscale(true).format("10^{%L}");
		//g.title("");
		//g.manual_settings("set pm3d map corners2color c1 \n set palette defined ( 0 \"white\", .01 \"blue\", 3 \"green\", 6 \"yellow\", 10 \"red\")");
		//g.plot3d(Plottable2DHistogram("",muonVsNeutrinoEnergy),"muonEnVsNeutrinoEn.eps");
		//g.manual_settings("");
		
		dataZenith.setUseContentScaling(false);
		/*std::cout << "dataZenith:\n" << dataZenith << std::endl;
		std::cout << "dataEnergy:\n" << dataEnergy << std::endl;
		std::cout << "convZenith:\n" << convZenith << std::endl;
		std::cout << "convEnergy:\n" << convEnergy << std::endl;
		std::cout << "promptZenith:\n" << promptZenith << std::endl;
		std::cout << "promptEnergy:\n" << promptEnergy << std::endl;
		std::cout << "astroZenith:\n" << astroZenith << std::endl;
		std::cout << "astroEnergy:\n" << astroEnergy << std::endl;*/
		
 		//compare CR flux slope
 		{
			g.title("");
 			auto softerWeight=DFWM(std::vector<double>{0.765801,0,0,.2,0,1,1,0});
 			auto defaultWeight=DFWM(std::vector<double>{0.765801,0,0,0,0,1,1,0});
 			auto harderWeight=DFWM(std::vector<double>{0.765801,0,0,-.2,0,1,1,0});
 			
 			Plottable1DHistogram softEnergy("Softer CR spectrum",Format::LEFT_STEPS,enAxis), 
 			                     normEnergy("Default CR spectrum",Format::LEFT_STEPS,enAxis), 
 			                     hardEnergy("Harder CR spectrum",Format::LEFT_STEPS,enAxis);
 			softEnergy.setUseContentScaling(false);
 			normEnergy.setUseContentScaling(false);
 			hardEnergy.setUseContentScaling(false);
 			softEnergy.format().line_width(2);
 			normEnergy.format().line_width(2);
 			hardEnergy.format().line_width(2);
 			for(const Event& e : mainSimulation){
 				softEnergy.add(e.energy,amount(softerWeight(e)));
 				normEnergy.add(e.energy,amount(defaultWeight(e)));
 				hardEnergy.add(e.energy,amount(harderWeight(e)));
 			}
 			
 			g.x_axis().label("Energy Proxy (arb. units)").logscale(true).format("10^{%L}").min(2e2);
 			g.y_axis().label("Events").logscale(true).format("10^{%L}").auto_max().min(.1);
 			
 			g.plot(softEnergy,normEnergy,hardEnergy,plot_path+"CR_slope_var.eps");
 		}
 		
 		//compare DOM eff
 		{
 			auto lowWeight=DFWM(std::vector<double>{0.765801,0,0,0,.0,1,1,0});
 			auto defaultWeight=DFWM(std::vector<double>{0.765801,0,0,0,0,1,1,0});
 			auto highWeight=DFWM(std::vector<double>{0.765801,0,0,0,0.2,1,1,0});
 			
 			Plottable1DHistogram lowEnergy("Lower Detector Efficiency",Format::LEFT_STEPS,enAxis), 
 			                     normEnergy("Default Detector Efficiency",Format::LEFT_STEPS,enAxis), 
 			                     highEnergy("Higher Detector Efficiency",Format::LEFT_STEPS,enAxis);
 			lowEnergy.setUseContentScaling(false);
 			normEnergy.setUseContentScaling(false);
 			highEnergy.setUseContentScaling(false);
 			lowEnergy.format().line_width(2);
 			normEnergy.format().line_width(2);
 			highEnergy.format().line_width(2);
 			for(const Event& e : mainSimulation){
 				lowEnergy.add(0.9*e.energy,amount(lowWeight(e)));
 				normEnergy.add(e.energy,amount(defaultWeight(e)));
 				highEnergy.add(1.1*e.energy,amount(highWeight(e)));
 			}
 			
 			g.x_axis().label("Energy Proxy (arb. units)").logscale(true).format("10^{%L}").min(2e2);
 			g.y_axis().label("Events").logscale(true).format("10^{%L}").auto_max().min(.1);
 			
 			g.plot(lowEnergy,normEnergy,highEnergy,plot_path+"DOM_eff_var.eps");
 		}
		
		//compare pi/K
		{
			auto lowWeight=DFWM(std::vector<double>{0.765801,0,0,0,0,.5,1,0});
			auto defaultWeight=DFWM(std::vector<double>{0.765801,0,0,0,0,1,1,0});
			auto highWeight=DFWM(std::vector<double>{0.765801,0,0,0,0,1.5,1,0});
			
			Plottable1DHistogram lowEnergy("Lower Kaon Contribution",Format::LEFT_STEPS,enAxis),
			normEnergy("Default Kaon Contribution",Format::LEFT_STEPS,enAxis),
			highEnergy("Higher Kaon Contribution",Format::LEFT_STEPS,enAxis);
			lowEnergy.setUseContentScaling(false);
			normEnergy.setUseContentScaling(false);
			highEnergy.setUseContentScaling(false);
			lowEnergy.format().line_width(2);
			normEnergy.format().line_width(2);
			highEnergy.format().line_width(2);
			for(const Event& e : mainSimulation){
				lowEnergy.add(0.9*e.energy,amount(lowWeight(e)));
				normEnergy.add(e.energy,amount(defaultWeight(e)));
				highEnergy.add(1.1*e.energy,amount(highWeight(e)));
			}
			
			g.x_axis().label("Energy Proxy (arb. units)").logscale(true).format("10^{%L}").min(2e2);
			g.y_axis().label("Events").logscale(true).format("10^{%L}").auto_max().min(.1);
			
			g.plot(lowEnergy,normEnergy,highEnergy,plot_path+"Kaon_var.eps");
		}
	}
	
	if(!outputFile.empty())
		std::cout.rdbuf(stdoutBackup);
	
	return(0);
}
