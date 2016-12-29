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

namespace DOMEff1{
	using DECFq=DOMEffCorrectionFactory::quadratic;
	
	//IC79, Photonics
	/*
	DOMEffCorrectionFactory convDOMEffParams{DECFq{-6.9575,20.7311,-10.1683},//A
									   DECFq{-188.55,1063.3,-98.888},//Ec
									   DECFq{5e-07,-0.00065081,0.00136781},//a
									   DECFq{0.16429,-0.372024,0.138571},//b
									   DECFq{0.1519,-0.26381,1.0547},//c
									   DECFq{-4.53903,10.4832,-6.00185},//d
									   DECFq{-0.2315,0.37835,-5.11185}//f
	};

	DOMEffCorrectionFactory promptDOMEffParams{DECFq{-0.0021655,0.0133636,-0.00327631},//A
									   DECFq{559.3,-462.54,650.783},//Ec
									   DECFq{-0.00203415,0.00196067,0.000808943},//a
									   DECFq{-0.01763,0.023565,-0.0685176},//b
									   DECFq{0.06695,-0.104605,0.974406},//c
									   DECFq{0.58235,-0.493365,-1.14509},//d
									   DECFq{-0.158,2.1909,-7.23742}//f
	};

	DOMEffCorrectionFactory astroDOMEffParams{DECFq{-0.0003025,0.00174725,-0.00014423},//A
									   DECFq{713.95,-840.035,922.137},//Ec
									   DECFq{-0.0022845,0.00995945,-0.00423254},//a
									   DECFq{0.027675,-0.0838325,-0.0316132},//b
									   DECFq{-0.15255,0.041695,0.89142},//c
									   DECFq{0.4403,0.71739,-1.89955},//d
									   DECFq{2.393,-4.4035,-3.08024}//f
	};
	*/

	//IC79, PPC
	DOMEffCorrectionFactory convDOMEffParams{DECFq{-1.5095,9.21065,-5.11306},//A
											 DECFq{-2933.65,6204.47,-2556.21},//Ec
											 DECFq{0.253088,-0.470608,0.220566},//a
											 DECFq{0.364265,-1.32139,0.834971},//b
											 DECFq{-12.4538,27.7233,-13.9344},//c
											 DECFq{24.6677,-52.6563,26.9565},//d
											 DECFq{32.198,-60.7569,23.7681}//f
	};
	DOMEffCorrectionFactory promptDOMEffParams{DECFq{0.0021485,0.00588975,-0.00186549},//A
											   DECFq{-1256.6,2660.29,-699.703},//Ec
											   DECFq{-0.00181225,0.00440502,-0.00177189},//a
											   DECFq{0.025705,-0.0638095,-0.0231598},//b
											   DECFq{0.11305,-0.187285,1.00776},//c
											   DECFq{0.4385,-0.62355,-0.94293},//d
											   DECFq{42.6655,-84.0792,36.4446}//f
	};
	DOMEffCorrectionFactory astroDOMEffParams{DECFq{0.0033839,-0.00543089,0.00308648},//A
											  DECFq{868.4,-1529.67,1414.94},//Ec
											  DECFq{-0.003235,0.0092993,0.00108512},//a
											  DECFq{1.15715,-2.47491,1.16225},//b
											  DECFq{-11.8548,24.9322,-11.9012},//c
											  DECFq{33.2671,-68.708,33.9192},//d
											  DECFq{48.4945,-95.7944,42.4776}//f
	};
}

namespace DOMEff2{
	using DECFq=DOMEffCorrectionFactory::quadratic;
	
	DOMEffCorrectionFactory convDOMEffParams2010{
		DECFq::fit(1.,1.59248,1.1,1.95364,1.21,2.32851),
		DECFq::fit(1.,685.621,1.1,728.915,1.21,781.983),
		DECFq::fit(1.,0.0017509,1.1,0.00142394,1.21,0.00103492),
		DECFq::fit(1.,-0.0615943,1.1,-0.0618721,1.21,-0.0633076),
		DECFq::fit(1.,0.892393,1.1,0.894665,1.21,0.896102),
		DECFq::fit(1.,-0.172292,1.1,-0.09683,1.21,0.00259697),
		DECFq::fit(1.,-4.906,1.1,-5.06695,1.21,-5.15766),
		DECFq::fit(1.,14.5788,1.1,12.5419,1.21,10.7938),
		DECFq::fit(1.,2.78094,1.1,2.80578,1.21,2.83577),
		DECFq::fit(1.,0.12,1.1,0.12,1.21,0.12),
		DECFq::fit(1.,-20.0357,1.1,-18.247,1.21,-15.8098),
		DECFq::fit(1.,2.83006,1.1,2.86105,1.21,2.89076),
		DECFq::fit(1.,0.15,1.1,0.15,1.21,0.15),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};
	DOMEffCorrectionFactory promptDOMEffParams2010{
		DECFq::fit(1.,0.00537151,1.1,0.00537151,1.21,0.0065974),
		DECFq::fit(1.,684.795,1.1,721.209,1.21,756.223),
		DECFq::fit(1.,0.00090789,1.1,0.00093,1.21,0.000960934),
		DECFq::fit(1.,-0.0594536,1.1,-0.0608898,1.21,-0.0621554),
		DECFq::fit(1.,0.929986,1.1,0.9345,1.21,0.936374),
		DECFq::fit(1.,-1.16802,1.1,-1.14,1.21,-1.11857),
		DECFq::fit(1.,-4.95612,1.1,-5.09910,1.21,-5.26351),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};
	DOMEffCorrectionFactory astroDOMEffParams2010{
		DECFq::fit(1.,0.000434767,1.1,0.000456363,1.21,0.000472558),
		DECFq::fit(1.,712.452,1.1,744.388,1.21,766.749),
		DECFq::fit(1.,0.0033222,1.1,0.00509226,1.21,0.00707188),
		DECFq::fit(1.,-0.0979135,1.1,-0.132507,1.21,-0.168295),
		DECFq::fit(1.,0.935332,1.1,1.14243,1.21,1.34662),
		DECFq::fit(1.,-1.30503,1.1,-1.68485,1.21,-2.05482),
		DECFq::fit(1.,-4.85442,1.1,-4.92682,1.21,-5.35705),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};
	
	/*DOMEffCorrectionFactory convDOMEffParams2010{
		DECFq::fit(.9,1.95383,1.0,2.58809,1.1,3.19216),
		DECFq::fit(.9,651.565,1.0,714.619,1.1,719),
		DECFq::fit(.9,0.00202053,1.0,0.00304648,1.1,0.00913419),
		DECFq::fit(.9,-0.0592227,1.0,-0.122151,1.1,-0.177794),
		DECFq::fit(.9,0.928953,1.0,1.33505,1.1,1.49207),
		DECFq::fit(.9,-0.453255,1.0,-1.03201,1.1,-1.11741),
		DECFq::fit(.9,-4.83274,1.0,-4.79081,1.1,-4.10492),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};
	DOMEffCorrectionFactory promptDOMEffParams2010{
		DECFq::fit(.9,0.00517557,1.0,0.00617276,1.1,0.00721292),
		DECFq::fit(.9,676.712,1.0,703.987,1.1,706.13),
		DECFq::fit(.9,0.000724707,1.0,0.000820882,1.1,0.000880812),
		DECFq::fit(.9,-0.0597673,1.0,-0.0612643,1.1,-0.0622472),
		DECFq::fit(.9,0.930776,1.0,0.933527,1.1,0.938539),
		DECFq::fit(.9,-1.14894,1.0,-1.12798,1.1,-1.09825),
		DECFq::fit(.9,-4.66764,1.0,-4.96912,1.1,-4.41729),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};
	DOMEffCorrectionFactory astroDOMEffParams2010{
		DECFq::fit(.9,0.000939638,1.0,0.00103949,1.1,0.00120702),
		DECFq::fit(.9,741.643,1.0,753.672,1.1,783.069),
		DECFq::fit(.9,0.00683414,1.0,0.00714942,1.1,0.0074),
		DECFq::fit(.9,-0.127877,1.0,-0.15551,1.1,-0.16),
		DECFq::fit(.9,0.935363,1.0,1.17616,1.1,1.17986),
		DECFq::fit(.9,-0.971678,1.0,-1.52173,1.1,-1.40644),
		DECFq::fit(.9,-4.45688,1.0,-4.82237,1.1,-4.21797),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};*/
	
	DOMEffCorrectionFactory convDOMEffParams2011{
		DECFq::fit(1.,2.35715,1.1,2.62321,1.21,3.1986),
		DECFq::fit(1.,666.067,1.1,698.4,1.21,742.072),
		DECFq::fit(1.,0.00124305,1.1,0.00114632,1.21,0.0009),
		DECFq::fit(1.,-0.0625227,1.1,-0.0633522,1.21,-0.064),
		DECFq::fit(1.,0.900783,1.1,0.901122,1.21,0.91),
		DECFq::fit(1.,-0.0671324,1.1,-0.031074,1.21,0),
		DECFq::fit(1.,-4.53948,1.1,-4.82229,1.21,-5.04689),
		DECFq::fit(1.,20.9548,1.1,13.2515,1.21,12.377),
		DECFq::fit(1.,2.78007,1.1,2.78969,1.21,2.79732),
		DECFq::fit(1.,0.12,1.1,0.12,1.21,0.12),
		DECFq::fit(1.,-27.9288,1.1,-21.0023,1.21,-24.8666),
		DECFq::fit(1.,2.82664,1.1,2.83886,1.21,2.84433),
		DECFq::fit(1.,0.15,1.1,0.15,1.21,0.15),
		DECFq::fit(1.,3.7669,1.1,11.0338,1.21,23.3389),
		DECFq::fit(1.,3.22336,1.1,3.08792,1.21,3.06697),
		DECFq::fit(1.,0.27,1.1,0.27,1.21,0.27),
	};
	DOMEffCorrectionFactory promptDOMEffParams2011{
		DECFq::fit(1.,0.007228,1.1,0.0074792,1.21,0.009),
		DECFq::fit(1.,653.214,1.1,684.376,1.21,727.92),
		DECFq::fit(1.,0.000488451,1.1,0.000661715,1.21,0.00088),
		DECFq::fit(1.,-0.05966,1.1,-0.0614057,1.21,-0.0635685),
		DECFq::fit(1.,0.936319,1.1,0.937381,1.21,0.938789),
		DECFq::fit(1.,-1.10583,1.1,-1.08688,1.21,-1.0394),
		DECFq::fit(1.,-4.73473,1.1,-4.91103,1.21,-5.14629),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};
	DOMEffCorrectionFactory astroDOMEffParams2011{
		DECFq::fit(1.,0.000570719,1.1,0.000574135,1.21,0.0006),
		DECFq::fit(1.,678.342,1.1,719.097,1.21,742.504),
		DECFq::fit(1.,0.00359341,1.1,0.00359341,1.21,0.00359341),
		DECFq::fit(1.,-0.0995974,1.1,-0.0995974,1.21,-0.0995974),
		DECFq::fit(1.,0.90998,1.1,0.90998,1.21,0.90998),
		DECFq::fit(1.,-1.13589,1.1,-1.13589,1.21,-1.13589),
		DECFq::fit(1.,-4.72396,1.1,-4.75125,1.21,-5.00637),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};
	
	/*DOMEffCorrectionFactory convDOMEffParams2011{
		DECFq::fit(.9,1.95383,1.0,2.58809,1.1,3.19216),
		DECFq::fit(.9,651.565,1.0,714.619,1.1,719),
		DECFq::fit(.9,0.00202053,1.0,0.00304648,1.1,0.00913419),
		DECFq::fit(.9,-0.0592227,1.0,-0.122151,1.1,-0.177794),
		DECFq::fit(.9,0.928953,1.0,1.33505,1.1,1.49207),
		DECFq::fit(.9,-0.453255,1.0,-1.03201,1.1,-1.11741),
		DECFq::fit(.9,-4.83274,1.0,-4.79081,1.1,-4.10492),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};
	DOMEffCorrectionFactory promptDOMEffParams2011{
		DECFq::fit(.9,0.00517557,1.0,0.00617276,1.1,0.00721292),
		DECFq::fit(.9,676.712,1.0,703.987,1.1,706.13),
		DECFq::fit(.9,0.000724707,1.0,0.000820882,1.1,0.000880812),
		DECFq::fit(.9,-0.0597673,1.0,-0.0612643,1.1,-0.0622472),
		DECFq::fit(.9,0.930776,1.0,0.933527,1.1,0.938539),
		DECFq::fit(.9,-1.14894,1.0,-1.12798,1.1,-1.09825),
		DECFq::fit(.9,-4.66764,1.0,-4.96912,1.1,-4.41729),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};
	DOMEffCorrectionFactory astroDOMEffParams2011{
		DECFq::fit(.9,0.000939638,1.0,0.00103949,1.1,0.00120702),
		DECFq::fit(.9,741.643,1.0,753.672,1.1,783.069),
		DECFq::fit(.9,0.00683414,1.0,0.00714942,1.1,0.0074),
		DECFq::fit(.9,-0.127877,1.0,-0.15551,1.1,-0.16),
		DECFq::fit(.9,0.935363,1.0,1.17616,1.1,1.17986),
		DECFq::fit(.9,-0.971678,1.0,-1.52173,1.1,-1.40644),
		DECFq::fit(.9,-4.45688,1.0,-4.82237,1.1,-4.21797),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,0,1.1,0,1.21,0),
		DECFq::fit(1.,1,1.1,1,1.21,1),
	};*/
}

namespace DOMEff3{
	//std::unique_ptr<Splinetable> domEffConv2010(new Splinetable("dom_eff_fits/conv_IC79.fits"));
	//std::unique_ptr<Splinetable> domEffConv2011(new Splinetable("/home/carguelles/work/TheSterileSearch/LikelihoodFit/dom_eff_fits/conv_IC86.fits"));
	//std::unique_ptr<Splinetable> domEffConvPion2011(new Splinetable("/home/carguelles/work/TheSterileSearch/LikelihoodFit/dom_eff_fits/conv_pion_IC86.fits"));
	//std::unique_ptr<Splinetable> domEffConvKaon2011(new Splinetable("/home/carguelles/work/TheSterileSearch/LikelihoodFit/dom_eff_fits/conv_kaon_IC86.fits"));
	//std::unique_ptr<Splinetable> domEffPrompt2011(new Splinetable("/home/carguelles/work/TheSterileSearch/LikelihoodFit/dom_eff_fits/prompt_IC86.fits"));

  // new system
	std::unique_ptr<Splinetable> domEffConv2011;
	std::unique_ptr<Splinetable> domEffConvPion2011;
	std::unique_ptr<Splinetable> domEffConvKaon2011;
	std::unique_ptr<Splinetable> domEffPrompt2011;

  /*
	std::unique_ptr<Splinetable> domEffPrompt2010(new Splinetable("dom_eff_fits/prompt_IC79.fits"));
	std::unique_ptr<Splinetable> domEffAstro2010(new Splinetable("dom_eff_fits/astro_IC79.fits"));
	std::unique_ptr<Splinetable> domEffAstro2011(new Splinetable("dom_eff_fits/astro_IC86.fits"));
  */
}
