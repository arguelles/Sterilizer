#include <iterator>
#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <chrono>
#include <queue>

#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/optional.hpp>
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

#include "timer.h"
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

// Constructor
Sterilizer::Sterilizer(DataPaths dataPaths, SteeringParams steeringParams, SterileNeutrinoParameters snp){

  // Set the parameter sets of this object                               
  steeringParams_=steeringParams;
  dataPaths_=dataPaths;
  sterileNuParams=snp;

  // Set up the RNG                                                      
  SetRandomNumberGeneratorSeed(steeringParams.rngSeed);

  // Load the data
  if(readCompact_){
    LoadCompactData(dataPaths_.compact_file_path);
    LoadCompactMC(dataPaths_.compact_file_path);
  } else {
    LoadData(dataPaths_.data_path);
    LoadMC(dataPaths_.mc_path);
  }
  LoadFluxes(dataPaths_.flux_splines_path,snp);
}



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


struct fitResult{
	std::vector<double> params;
	double likelihood;
	unsigned int nEval, nGrad;
	bool succeeded;
};


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
	minimizer.addParameter(seed[6],.001,0.0,2.0);
	minimizer.addParameter(seed[7],.001,-1.0,1.0);
	
	for(auto idx : indicesToFix)
		minimizer.fixParameter(idx);
	
	minimizer.setChangeTolerance(1e-5);
	minimizer.setHistorySize(20);
	
	fitResult result;
	result.succeeded=minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));
	result.likelihood=minimizer.minimumValue();
	result.params=minimizer.minimumPosition();
	result.nEval=minimizer.numberOfEvaluations();
	result.nGrad=minimizer.numberOfEvaluations(); //gradient is always eval'ed with function
		
	return(result);
}



template<typename ProblemType>
double saturatedPoisson(const ProblemType& prob, const std::vector<double>& params){
	auto spProb=prob.makeAlternateLikelihood(likelihood::saturatedPoissonLikelihood());
	return(-spProb.evaluateLikelihood(params));
}


//precompute weights as much as possible
/////////////////////////////////////////////////// VERY IMPORTANT // MUY IMPORTANTE ///////////////////////////////////////////////////
template<typename ContainerType, typename WeighterType>
void initializeSimulationWeights(ContainerType& simulation, const WeighterType& convPionWeighter, const WeighterType& convKaonWeighter, const WeighterType& promptWeighter, const OversizeWeighter& osw){
	using iterator=typename ContainerType::iterator;
  //std::cout << simulation.size() << std::endl;
	auto cache=[&](iterator it, iterator end){
    //std::cout << "alala" << std::endl;
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
			e.cachedPromptWeight=0.;
      //e.cachedPromptWeight=promptWeighter(lw_e)*e.cachedLivetime;
      //std::cout << osweight << " " << e.cachedConvPionWeight << " " << e.cachedConvKaonWeight << " " << e.cachedPromptWeight << std::endl;
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

  std::string model_name = "HondaGaisser";  
  std::string delta_model_name = "HondaGaisser";


  std::pair<unsigned int, double> fixedParams;
  std::vector<double> fitSeed{1.02,0,0,0.05,.0985,1.1,1,0};
  std::vector<double> dataChallenge_nuisance_parameters{1,0,0,0,.1,1,1,0};
  std::vector<double> data{1,0,0,0,.1,1,1,0};
  std::vector<double> existingBest_(existingBest);
  std::vector<double> fitSeed_(fitSeed);
  std::vector<double> dataChallenge_nuisance_parameters_(dataChallenge_nuisance_parameters);
  unsigned int number_of_data_challenges = 1;
  std::vector<double> dataChallenge_seeds {1234};
  std::vector<double> dataChallenge_seeds_(dataChallenge_seeds);
  bool use_datachallenge_histogram=false;

	bool doDataChallenge=false;
	int yearsIC86=1;

  std::string dcsimulation_to_load = "nufsgen_mie_0_99";
  std::string xs_model_name = "";
  std::string data_challenge_histogram_filepath = "";
  std::string modelId = "";

	OptionParser op;
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
	op.addOption("dataChallenge",[&](){doDataChallenge=true;},"Run fit trials on MC realizations.");
	op.addOption("UseFactorization",[&](){use_factorization_technique=true;},"Will use factorized oscillations and propagation.");
	op.addOption("save_flux",[&](){save_flux=true;},"Save flux..");
	op.addOption("yearsIC86",yearsIC86,"Number of years of IC86 equivalent livetime.");
  op.addOption("th24_null",th24_null,"th24 sterile neutrino mixing paramater null hypothesis [rad].");
  op.addOption("th14",th14,"th14 sterile neutrino mixing paramater [rad].");
  op.addOption("th24",th24,"th24 sterile neutrino mixing paramater [rad].");
  op.addOption("th34",th34,"th34 sterile neutrino mixing paramater [rad].");
  op.addOption("del14",del14,"del14 sterile neutrino mixing paramater [rad].");
  op.addOption("del24",del24,"del24 sterile neutrino mixing paramater [rad].");
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
  op.addOption("flux_splines_path",flux_splines_path,"Path to the propaged flux splines");
  op.addOption("output_path",output_path, "String that specifies the output data path");
  op.addOption("oversize_function_main",oversize_function_main,"Main oversize correction function to use");
  op.addOption("oversize_function_dc",oversize_function_dc,"Data challenge oversize correction function to use");
  op.addOption("number_of_data_challenges",number_of_data_challenges,"Integer specifying the number of data challenges to perform.");
  op.addOption("dataChallenge_seeds",dataChallenge_seeds_,"Seeds used for each the data challenges.");
  op.addOption("modelId",modelId,"Sterile model id based on scan_values file.");
  op.addOption("use_datachallenge_histogram",[&](){use_datachallenge_histogram = true;},"If use_datachallenge_histogram is set to true then the DC MC would not be loaded and the histogram will be loaded instead.");
  op.addOption("data_challenge_histogram_filepath",data_challenge_histogram_filepath,"Path to the histogram use for data challenge when use_datachallenge_histogram is set to true.");
  op.addOption("save_datachallenge_histogram",[&](){save_dc_sample = true;},"If save_dc_sample is set to true then the DC MC will be saved with filename data_challenge_histogram_filepath.");
  op.addOption("do_asimov_sensitivity",[&](){do_asimov_sensitivity = true;},"if do_asimov_sensitivity is set to true then asimov sensitivity will be calculated.");

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
      if( e.energy > 4.e2 and e.energy < 2.e4 ){
      //{
        //dump_data << cos(e.zenith) << ' ' << e.azimuth << ' ' << e.energy << ' ';
        dump_data << e.energy << ' ' << e.zenith << ' ';
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


/////////////////////////////////////////////////// VERY IMPORTANT // MUY IMPORTANTE ///////////////////////////////////////////////////
// HERE WE CONSTRUCT THE LW WEIGHTERS FOR THE OSCILLATION HYPOTHESIS AND WEIGH THE MC //
// AQUI CONSTRUIMOS LOS PESADORES DE LW PARA UNA HIPOTESIS DE OSCILACION Y PESAMOS EL MC //
	t.start("WeightingMC");
  // flux weighters // pesadores del flujo
  if(!quiet)
    std::cout << "Begin Loading nuSQuIDS objects." << std::endl;

  std::shared_ptr<LW::Flux> flux_kaon,flux_pion;
  std::string sterile_neutrino_model_identifier = modelId+"_"+std::to_string(dm41sq)+"_"+std::to_string(th14)+"_"+std::to_string(th24)+"_"+std::to_string(th34)+"_"+std::to_string(del14)+"_"+std::to_string(del24);

  if(use_factorization_technique){
    std::string oscillation_model = "noint_atmospheric_"+sterile_neutrino_model_identifier+".hdf5";
    std::string atmospheric_model_pion = model_name + "_pion";
    std::string atmospheric_model_kaon = model_name + "_kaon";

    flux_pion = std::make_shared<LW::FactorizedSQUIDSFlux>(squids_files_path+oscillation_model,
                                                           flux_splines_path+atmospheric_model_pion+"_neutrino_spline.fits",
                                                           flux_splines_path+atmospheric_model_pion+"_antineutrino_spline.fits");
    flux_kaon = std::make_shared<LW::FactorizedSQUIDSFlux>(squids_files_path+oscillation_model,
                                                           flux_splines_path+atmospheric_model_kaon+"_neutrino_spline.fits",
                                                           flux_splines_path+atmospheric_model_kaon+"_antineutrino_spline.fits");
  } else{
      std::string flux_pion_filename = "pion_atmospheric_"+sterile_neutrino_model_identifier;
      std::string flux_kaon_filename = "kaon_atmospheric_"+sterile_neutrino_model_identifier;
      if(model_name != ""){
        flux_pion_filename+="_"+model_name;
        flux_kaon_filename+="_"+model_name;
      }
      flux_kaon = std::make_shared<LW::SQUIDSFlux>(squids_files_path + flux_kaon_filename + ".hdf5");
      flux_pion = std::make_shared<LW::SQUIDSFlux>(squids_files_path + flux_pion_filename + ".hdf5");
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
  //std::cout << "begin: initialize simulation weights" << std::endl;
	initializeSimulationWeights(mainSimulation,PionFluxWeighter,KaonFluxWeighter,PromptFluxWeighter,osw_main);
  //std::cout << "end: initialize simulation weights" << std::endl;

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
    auto weighter=DFWM(dataChallenge_nuisance_parameters);
    dump_mc << std::scientific;
    dump_mc.precision(6);
    dump_mc.fill('0');
    for(const Event& e : mainSimulation){
      if( e.energy > 4.e2 and e.energy < 2.e4 ){
      //{
        auto w=weighter(e);
        dump_mc << static_cast<int>(e.primaryType) << ' ';
        dump_mc << e.energy << ' ' << e.zenith << ' ';
        dump_mc << e.injectedEnergy << ' ' << e.intX << ' ' << e.intY << ' ';
        dump_mc << e.injectedMuonEnergy << ' ' << e.injectedMuonZenith << ' ';
        dump_mc << e.inelasticityProbability << ' ' << e.totalColumnDepth << ' ';
     //   dump_mc << e.cachedConvKaonWeight << std::endl;
        dump_mc << w << std::endl;
      }
    }
    dump_mc.close();
    std::cout << "End Dump MC info." << std::endl;
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
	dataHist.getAxis(0)->setLowerLimit(minFitEnergy);
	dataHist.getAxis(0)->setUpperLimit(maxFitEnergy);
  dataHist.getAxis(1)->setLowerLimit(minCosth);
  dataHist.getAxis(1)->setUpperLimit(maxCosth);
	
  if(!quiet)
    std::cout << "Fill data histogram" << std::endl;
	bin(burnsample,dataHist,binner);
	
	decltype(mainSimulation) simSubset1, simSubset2;

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
                                            {"HondaGaisser",8./7.},
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
          flux_kaon_null = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "kaon_atmospheric_"+modelId+"_"+std::to_string(dm41sq_null)+"_"+std::to_string(th24_null)+".hdf5");
          flux_pion_null = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "pion_atmospheric_"+modelId+"_"+std::to_string(dm41sq_null)+"_"+std::to_string(th24_null)+".hdf5");
        } else {
          flux_kaon_null = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "kaon_atmospheric_"+modelId+"_"+std::to_string(dm41sq_null)+"_"+std::to_string(th24_null)+"_"+model_name+".hdf5");
          flux_pion_null = std::make_shared<LW::SQUIDSFlux>(squids_files_path + "pion_atmospheric_"+modelId+"_"+std::to_string(dm41sq_null)+"_"+std::to_string(th24_null)+"_"+model_name+".hdf5");
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

        // esto necesita ser revisado. Hay un problema aqui.
        if(do_asimov_sensitivity){
          if(!quiet)
            std::cout << "Doing Asimov-sensitivity" << std::endl;
          if(is_same_simulation)
            bin(mainSimulation,sampleHist,binner);
          else
            bin(DCSimulation,sampleHist,binner);
        }
        else {
          std::vector<Event> sample;
          if(is_same_simulation)
            sample=likelihood::generateSample(weights,mainSimulation,expected,rng);
          else
            sample=likelihood::generateSample(weights,DCSimulation,expected,rng);
          if(!quiet)
            std::cout << " sampled " << sample.size() << " events" << std::endl;
          bin(sample,sampleHist,binner);
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
      fr = doFitLBFGSB(prob,seed,fixedIndices);
      
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


    t.stop();
    if(!quiet)
      std::cout << t << std::endl;
    return(0);


const std::string sterile_params_str = "_"+modelId+"_"+std::to_string(dm41sq)+"_"+std::to_string(th14)+"_"+std::to_string(th24)+"_"+std::to_string(th34)+"_"+std::to_string(del14)+"_"+std::to_string(del24)+".hdf5";


// captain - full stop
 exit(0);

}



	
//Set the RNG seed
void Sterilizer::SetRandomNumberGeneratorSeed(unsigned int seed)
{
  steeringParams_.rngSeed=seed;
  rng_.seed(seed);
  if(!steeringParams_.quiet) std::cout<<"setting RNG seed to " << seed<<std::endl;
}
 

// Check that the directories where files are mean to be exist               
 bool Sterilizer::CheckDataPaths(DataPaths dp)
 {
   CheckDataPath(dp.compact_file_path);
   CheckDataPath(dp.squids_files_path);
   CheckDataPath(dp.prompt_squids_files_path);
   CheckDataPath(dp.xs_spline_path);
   CheckDataPath(dp.data_path);
   CheckDataPath(dp.mc_path);
   CheckDataPath(dp.oversize_function_path);
   CheckDataPath(dp.domeff_spline_path);
   CheckDataPath(dp.flux_splines_path);
 }

 // Check a directory exists and throw a relevant error otherwise.            
 bool Sterilizer::CheckDataPath(std::string p)
 {
   struct stat info;
   if(p!="")
     {
       if( stat( p, &info ) != 0 )
	 {
	   throw std::runtime_error("cannot access "+ p);
	 }
       else if( !(info.st_mode & S_IFDIR) )
	 {
	   throw std::runtime_error("is not a directory: " +p);
	 }
     }
   else
     std::cout<<"Warning, there are unset paths in DataPaths. Check you want \
this."<<std::endl;
 }
