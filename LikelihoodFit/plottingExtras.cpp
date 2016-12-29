#include "plottingExtras.h"

interval<unsigned int> fcPoissonAcceptanceInterval(double mean, double confidence, double offset){
	assert(confidence>0.0 && confidence<1.0);
	assert(offset>=0.0);
	
	auto poisson=[offset](unsigned int n, double mu)->double{
		return(exp((n*log(mu+offset)-lgamma(n+1)-(mu+offset))));
	};
	auto fcR=[offset](unsigned int n, double mu)->double{
		double muBest=std::max(0.0,n-offset);
		return(exp(n*(log(mu+offset)-log(muBest+offset))+muBest-mu));
	};
	
	//start from the n with the highest R, which will be one of the two closest to mean+offset
	interval<unsigned int> result;
	if(fcR(ceil(mean+offset),mean)>=fcR(floor(mean+offset),mean))
		result.min=result.max=ceil(mean+offset);
	else
		result.min=result.max=floor(mean+offset);
	//keep track of probability accumulated
	double pAcc=poisson(result.min,mean);
	
	//expand the interval outwards until enough probability has accumulated
	//keep track of the R values for the two values of n just outside the interval so that
	//we can always pick out the better one without having to repeat calculations of R
	unsigned int nextLower=std::max((int)0,(int)result.min-1);
	unsigned int nextHigher=result.max+1;
	//use -1 as an out-of-band value to indicate 'not yet computed'
	double nextLowerR=-1, nextHigherR=-1;
	while(pAcc<confidence){
		//make sure both candidates R values are computed
		if(nextLowerR==-1 && result.min>0)
			nextLowerR=fcR(nextLower,mean);
		if(nextHigherR==-1)
			nextHigherR=fcR(nextHigher,mean);
		
		//add the next best value of n into the interval
		if(nextHigherR>=nextLowerR){
			result.max=nextHigher;
			nextHigher++;
			nextHigherR=-1;
			pAcc+=poisson(result.max,mean);
		}
		else{
			result.min=nextLower;
			if(nextLower>0)
				nextLower--;
			nextLowerR=-1;
			pAcc+=poisson(result.min,mean);
		}
	}
	
	return(result);
}

interval<double> fcPoissonConfidenceInterval(unsigned int obs, double confidence, double offset, double tol){
	assert(offset>=0.0);
	
	//work outward from a value definitely in the interval
	double inner=obs-offset;
	interval<double> result;
	
	//find left side
	//bracketing is trivial
	double min=0.0, max=inner;
	//bisect
	while((max-min)>max*tol){
		double test=(min+max)/2;
		interval<unsigned int> accept=fcPoissonAcceptanceInterval(test, confidence, offset);
		(accept.contains(obs) ? max : min)=test;
	}
	result.min=(min+max)/2;
	//find right side
	//bracket
	min=inner;
	max=(obs?2*inner:1);
	while(fcPoissonAcceptanceInterval(max, confidence, offset).contains(obs))
		max*=2;
	//bisect
	while((max-min)>max*tol){
		double test=(min+max)/2;
		interval<unsigned int> accept=fcPoissonAcceptanceInterval(test, confidence, offset);
		(accept.contains(obs) ? min : max)=test;
	}
	result.max=(min+max)/2;
	
	return(result);
}

fcErrorValue operator*(double scale, const fcErrorValue& v){
	return(fcErrorValue(v)*=scale);
}

std::ostream& operator<<(std::ostream& os, const fcErrorValue& v){
	os << (double)v << '[' << v.errorMin() << ',' << v.errorMax() << ']';
	return(os);
}

sqErrorValue operator*(double scale, const sqErrorValue& v){
	return(sqErrorValue(v)*=scale);
}

std::ostream& operator<<(std::ostream& os, const sqErrorValue& v){
	os << (double)v << '[' << v.errorMin() << ',' << v.errorMax() << ']';
	return(os);
}