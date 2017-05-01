#include "SterileSearch.h"
#include "GenerationSpecifications.h"

namespace SterileSearch {

/*************************************************************************************************************
 * Constructor
 * **********************************************************************************************************/

Sterilizer::Sterilizer(DataPaths dataPaths, SteeringParams steeringParams, SterileNuParams snp):
  steeringParams_(steeringParams),dataPaths_(dataPaths),sterileNuParams_(snp){
  if(!steeringParams_.quiet) std::cout<<"Sterilizer constructor: checking paths" <<std::endl;
  CheckDataPaths(dataPaths_);

  if(!steeringParams_.quiet) std::cout<<"Loading Splines" <<std::endl;
  LoadDOMEfficiencySplines();

  if(steeringParams_.ReadCompact){
    if(!steeringParams_.quiet) std::cout<<"Loading compact data" <<std::endl;
    LoadCompact();
  } else {
    if(!steeringParams_.quiet) std::cout<<"Loading data" <<std::endl;
    LoadData();
    if(!steeringParams_.quiet) std::cout<<"Loading MC" <<std::endl;
    LoadMC();
  }
  if(!steeringParams_.quiet) std::cout<<"Loading  XS" <<std::endl;
  ConstructCrossSectionWeighter();
  if(!steeringParams_.quiet) std::cout<<"Loading Flux weighter" <<std::endl;
  ConstructFluxWeighter();
  if(!steeringParams_.quiet) std::cout<<"Loading MC weighter" <<std::endl;
  ConstructMonteCarloGenerationWeighter();
  if(!steeringParams_.quiet) std::cout<<"Loading Lepton weighter" <<std::endl;
  ConstructLeptonWeighter();
  if(!steeringParams_.quiet) std::cout<<"Loading Oversize weighter" <<std::endl;
  ConstructOversizeWeighter();
  if(!steeringParams_.quiet) std::cout<<"Weighting MC" <<std::endl;
  WeightMC();
  if(!steeringParams_.quiet) std::cout<<"Making data hist" <<std::endl;
  ConstructDataHistogram();
  if(!steeringParams_.quiet) std::cout<<"Making sim hist" <<std::endl;
  ConstructSimulationHistogram();
  if(!steeringParams_.quiet) std::cout<<"Construcing likelihood problem with default settings" <<std::endl;
  ConstructLikelihoodProblem(Priors(), Nuisance(),NuisanceFlag());
}

/*************************************************************************************************************
 * Implementation auxiliary functions
 * **********************************************************************************************************/

std::string GetSterileNeutrinoModelIdentifier(SterileNuParams snp){
  return std::to_string(snp.modelId)+"_"+std::to_string(snp.dm41sq)+"_"+std::to_string(snp.th14)+"_"+std::to_string(snp.th24)+"_"+std::to_string(snp.th34)+"_"+std::to_string(snp.del14)+"_"+std::to_string(snp.del24);
}

auto binner = [](HistType& h, const Event& e){
                h.add(e.energy,cos(e.zenith),e.year,amount(std::cref(e)));
};




/*************************************************************************************************************
 * Functions to read and write data
 * **********************************************************************************************************/

void Sterilizer::LoadData(){
  try{
    auto dataAction = [&](RecordID id, Event& e, int dataYear){
      if(e.check(false,Level::neutrino)){
        e.year=dataYear;
        e.cachedWeight=1.;
        sample_.push_back(e);
      }
    };
    auto ic86Action=[&](RecordID id, Event& e){ dataAction(id,e,2011); };
    if (steeringParams_.useBurnSample)
      readFile(CheckedFilePath(dataPaths_.data_path+"burnsample_ic86.h5"),ic86Action);
    else
      readFile(CheckedFilePath(dataPaths_.data_path+"IC86.h5"),ic86Action);    
  } catch(std::exception& ex){
    std::cerr << "Problem loading experimental data: " << ex.what() << std::endl;
  }
  if(!steeringParams_.quiet)
    std::cout << "Loaded " << sample_.size() << " experimental events" << std::endl;
  data_loaded_=true;
}


void Sterilizer::LoadMC(){

  if(not dom_efficiency_splines_constructed_)
    throw std::runtime_error("MC cannot be loaded until dom splines are loaded");
    std::map<unsigned int,double> livetime;
    if(!steeringParams_.useBurnSample)
      livetime=steeringParams_.fullLivetime;
    else
      livetime=steeringParams_.burnSampleLivetime;

    std::vector<std::string> simSetsToLoad;
    simSetsToLoad.push_back(steeringParams_.simToLoad.c_str());
    std::map<std::string,run> simInfo=GetSimInfo(dataPaths_.mc_path);
    
    struct domEffSetter{
      simpleEffRate<Event> convDOMEffRate;
      domEffSetter(double simulatedDOMEfficiency, std::shared_ptr<Splinetable> domEffConv):
	convDOMEffRate(domEffConv.get(),simulatedDOMEfficiency,&Event::cachedConvDOMEff) {}
      void setCache(Event& e) const{
	convDOMEffRate.setCache(e);
      }
    };

    try{
      auto simAction=[&](RecordID id, Event& e, int simYear, const domEffSetter& domEff){
	if(e.check(true,Level::neutrino) && e.energy>1){
	  e.year=simYear;
	  e.cachedLivetime=livetime.find(simYear)->second;
	  e.cachedConvPionWeight=0;
	  e.cachedConvKaonWeight=0;
	  e.cachedPromptWeight=0;
	  domEff.setCache(e);
	  if(e.primaryType==particleType::NuTau || e.primaryType==particleType::NuTauBar){
	    assert(e.cachedConvPionWeight==0.0);
	    assert(e.cachedConvKaonWeight==0.0);
	    assert(e.cachedPromptWeight==0.0);
	  }
	  mainSimulation_.push_back(e);
	}
      };

      for(auto simSet : simSetsToLoad){
	const auto& setInfo=simInfo.find(simSet)->second;
	unsigned int simYear=setInfo.details.year;
	//unsigned int yearindex=yearindices_[simYear];
	domEffSetter domEff(setInfo.unshadowedFraction,domEffConv_[simYear]);
	auto callback=[&,simYear](RecordID id, Event& e){ simAction(id,e,simYear,domEff); };
	auto path=CheckedFilePath(dataPaths_.mc_path+setInfo.filename);
	readFile(path,callback);
      }
    } catch(std::exception& ex) 
      {
	std::cerr << "Problem loading simulated data: " << ex.what() << std::endl;
      }
    if(!steeringParams_.quiet)
      std::cout << "Loaded " << mainSimulation_.size() << " events in main simulation set" << std::endl;
    simulation_loaded_=true;
}

void Sterilizer::LoadCompact(){
  std::map<std::string,run> simInfo=GetSimInfo(dataPaths_.mc_path);
  try{
    auto simulation_information=simInfo.find(steeringParams_.simToLoad);
    if(simulation_information==simInfo.end())
      throw std::runtime_error("Could not find " + steeringParams_.simToLoad + " in simulation list");

    std::cout<<"There is a segfault between this line"<<std::endl;
    std::cout<<std::string(dataPaths_.data_path+simulation_information->second.filename)<<std::endl;
    std::cout<<CheckedFilePath(std::string(dataPaths_.data_path+simulation_information->second.filename))<<std::endl;
    std::string original_data_file = CheckedFilePath(dataPaths_.data_path+"IC86.h5");
    std::cout<<" and this one"<<std::endl;

    std::string original_simulation_file = CheckedFilePath(dataPaths_.mc_path+(*simulation_information).second.filename);
    CheckedFilePath(dataPaths_.compact_file_path+"/"+steeringParams_.simToLoad+"_compact_data.dat");
    unsplatData(CheckedFilePath(dataPaths_.compact_file_path+"/"+steeringParams_.simToLoad+"_compact_data.dat"),
        getFileChecksum(original_data_file)+getFileChecksum(original_simulation_file),sample_,mainSimulation_);

    if(!steeringParams_.quiet){
      std::cout << "Loaded " << sample_.size() << " experimental events." << std::endl;
      std::cout << "Loaded " << mainSimulation_.size() << " events in main simulation set." << std::endl;
    }
  } catch(std::runtime_error& re){
    std::cerr << re.what() << std::endl;
    std::cerr << "Failed to load compact data" << std::endl;
  }
  data_loaded_=true;
  simulation_loaded_=true;
  std::cout<<"end of load compact"<<std::endl;
}

bool Sterilizer::WriteCompact() const {
  std::map<std::string,run> simInfo=GetSimInfo(dataPaths_.mc_path);
  try{
    auto simulation_information=simInfo.find(steeringParams_.simToLoad);
    if(simulation_information==simInfo.end())
      throw std::runtime_error("Could not find " + steeringParams_.simToLoad + " in simulation list");
    std::string original_data_file =dataPaths_.data_path+"IC86.h5";
    std::string original_simulation_file = dataPaths_.mc_path+(*simulation_information).second.filename;
    std::cout<<"Files " << original_data_file<<" " << original_simulation_file<<std::endl;
    splatData(dataPaths_.compact_file_path+"/"+steeringParams_.simToLoad+"_compact_data.dat",
	      getFileChecksum(original_data_file)+getFileChecksum(original_simulation_file),sample_,mainSimulation_);
    return true;
  } catch(std::runtime_error& re){
    std::cerr << re.what() << std::endl;
    std::cerr << "Failed to save compact data" << std::endl;
    return false;
  }
}

void Sterilizer::ClearData(){
  sample_.clear();
}

void Sterilizer::ClearSimulation(){
  mainSimulation_.clear();
}

/*************************************************************************************************************
 * Functions to load to load DOM efficiency splines
 * **********************************************************************************************************/

void Sterilizer::LoadDOMEfficiencySplines(){
  for(unsigned int year : steeringParams_.years){
    domEffConv_.insert({year,std::unique_ptr<Splinetable>(new Splinetable(CheckedFilePath(dataPaths_.domeff_spline_path+"/conv_IC"+std::to_string(year)+".fits")))});
    domEffPrompt_.insert({year,std::unique_ptr<Splinetable>(new Splinetable(CheckedFilePath(dataPaths_.domeff_spline_path+"/prompt_IC"+std::to_string(year)+".fits")))});
  }
  DFWM.SetSplines(domEffConv_,domEffPrompt_);
  dom_efficiency_splines_constructed_=true;
}

/*************************************************************************************************************
 * Functions to construct weighters
 * **********************************************************************************************************/

void Sterilizer::ConstructCrossSectionWeighter(){
  if(steeringParams_.xs_model_name=="")
    {  
      CheckedFilePath( dataPaths_.xs_spline_path + "/dsdxdy-numu-N-cc.fits");
      CheckedFilePath( dataPaths_.xs_spline_path + "/dsdxdy-numubar-N-cc.fits");
    }
  else
    {
      CheckedFilePath( dataPaths_.xs_spline_path + "/dsdxdy-numu-N-cc-"+steeringParams_.xs_model_name+".fits");
      CheckedFilePath( dataPaths_.xs_spline_path + "/dsdxdy-numubar-N-cc-"+steeringParams_.xs_model_name+".fits");
    }

  xsw_ = std::make_shared<LW::CrossSectionFromSpline>(static_cast<std::string>(dataPaths_.xs_spline_path),steeringParams_.xs_model_name);
  xs_weighter_constructed_=true;
}

void Sterilizer::ConstructFluxWeighter(){
  std::string sterile_neutrino_model_identifier = GetSterileNeutrinoModelIdentifier(sterileNuParams_);

  if(steeringParams_.useFactorization){
    std::string oscillation_model = "noint_atmospheric_"+sterile_neutrino_model_identifier+".hdf5";
    std::string atmospheric_model_pion = steeringParams_.modelName + "_pion";
    std::string atmospheric_model_kaon = steeringParams_.modelName + "_kaon";

    fluxPion_ = std::make_shared<LW::FactorizedSQUIDSFlux>(CheckedFilePath(dataPaths_.squids_files_path + oscillation_model),
                                                           CheckedFilePath(dataPaths_.flux_splines_path+atmospheric_model_pion+"_neutrino_spline.fits"),
                                                           CheckedFilePath(dataPaths_.flux_splines_path+atmospheric_model_pion+"_antineutrino_spline.fits"));
    fluxKaon_ = std::make_shared<LW::FactorizedSQUIDSFlux>(CheckedFilePath(dataPaths_.squids_files_path+oscillation_model),
                                                           CheckedFilePath(dataPaths_.flux_splines_path+atmospheric_model_kaon+"_neutrino_spline.fits"),
                                                           CheckedFilePath(dataPaths_.flux_splines_path+atmospheric_model_kaon+"_antineutrino_spline.fits"));
  } else {
      std::string flux_pion_filename = "pion_atmospheric_"+sterile_neutrino_model_identifier;
      std::string flux_kaon_filename = "kaon_atmospheric_"+sterile_neutrino_model_identifier;
      flux_pion_filename+="_"+steeringParams_.modelName;
      flux_kaon_filename+="_"+steeringParams_.modelName;

      if(dataPaths_.use_simple_filename)
	{
	  flux_pion_filename = "pion_"+std::to_string(sterileNuParams_.modelId);
	  flux_kaon_filename = "kaon_"+std::to_string(sterileNuParams_.modelId);
	}

      fluxKaon_ = std::make_shared<LW::SQUIDSFlux>(CheckedFilePath(dataPaths_.squids_files_path + flux_kaon_filename + ".hdf5"));
      fluxPion_ = std::make_shared<LW::SQUIDSFlux>(CheckedFilePath(dataPaths_.squids_files_path + flux_pion_filename + ".hdf5"));
  }
  std::string PromptPath;
  if(steeringParams_.onePromptFitsAll)
    PromptPath=CheckedFilePath(dataPaths_.prompt_squids_files_path + "prompt_atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5");
  else
    PromptPath=CheckedFilePath(dataPaths_.prompt_squids_files_path + "prompt_atmospheric_"+sterile_neutrino_model_identifier+".hdf5");
  fluxPrompt_ = std::make_shared<LW::SQUIDSFlux>(PromptPath);

  flux_weighter_constructed_=true;
}

void Sterilizer::ConstructMonteCarloGenerationWeighter(){
  std::map<std::string,run> simInfo=GetSimInfo(dataPaths_.mc_path);
  std::vector<std::string> simSetsToLoad;
  simSetsToLoad.push_back(steeringParams_.simToLoad);
  for( std::string sim_name : simSetsToLoad )
    mcw_.addGenerationSpectrum(simInfo.find(sim_name)->second.details);
  mc_generation_weighter_constructed_=true;
}

void Sterilizer::ConstructLeptonWeighter(){
  if(not mc_generation_weighter_constructed_)
    throw std::runtime_error("MonteCarlo generation weighter has to be constructed first.");
  if(not flux_weighter_constructed_)
    throw std::runtime_error("Flux weighter has to be constructed first.");
  if(not xs_weighter_constructed_)
    throw std::runtime_error("Cross section weighter has to be constructed first.");

  pionFluxWeighter_ = std::make_shared<LW::LeptonWeighter>(fluxPion_,xsw_,mcw_);
  kaonFluxWeighter_ = std::make_shared<LW::LeptonWeighter>(fluxKaon_,xsw_,mcw_);
  promptFluxWeighter_ = std::make_shared<LW::LeptonWeighter>(fluxPrompt_,xsw_,mcw_);
  lepton_weighter_constructed_=true;
}

void Sterilizer::ConstructOversizeWeighter(){
  osw_=OversizeWeighter(CheckedFilePath(dataPaths_.oversize_path+"/"+steeringParams_.oversizeFunction+".dat"));
  oversize_weighter_constructed_=true;
}

/*************************************************************************************************************
 * Functions to initialize the MC weights
 * **********************************************************************************************************/

void Sterilizer::WeightMC(){
  if(not lepton_weighter_constructed_)
    throw std::runtime_error("LeptonWeighter has to be constructed first.");
  if(not oversize_weighter_constructed_)
    throw std::runtime_error("OversizeWeighter has to be constructed first.");
  if(not simulation_loaded_)
    throw std::runtime_error("No simulation has been loaded. Cannot construct simulation histogram.");
  InitializeSimulationWeights();
  simulation_initialized_=true;
}

void Sterilizer::InitializeSimulationWeights(){
  if(not dom_efficiency_splines_constructed_)
    throw std::runtime_error("Simulation cannot be weighted until dom splines are loaded");
  using iterator=std::deque<Event>::iterator;
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
	  0, // for azimuth
      static_cast<particleType>(e.primaryType),
      static_cast<int>(e.year)};

      double osweight = osw_.EvaluateOversizeCorrection(e.energy, e.zenith);
      e.cachedConvPionWeight=(*pionFluxWeighter_)(lw_e)*e.cachedLivetime*osweight;
      e.cachedConvKaonWeight=(*kaonFluxWeighter_)(lw_e)*e.cachedLivetime*osweight;
      e.cachedPromptWeight=(*promptFluxWeighter_)(lw_e)*e.cachedLivetime*osweight;

    }
  };

  iterator it=mainSimulation_.begin(), end=mainSimulation_.end();
  cache(it,end);
  /*
  ThreadPool pool(steeringParams_.evalThreads);
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
  */
}

/*************************************************************************************************************
 * Functions to construct histograms
 * **********************************************************************************************************/

void Sterilizer::ConstructDataHistogram(){
  if(not data_loaded_)
    throw std::runtime_error("No data has been loaded. Cannot construct data histogram.");

  dataHist_ = HistType(LogarithmicAxis(steeringParams_.logEbinEdge, steeringParams_.logEbinWidth), LinearAxis(steeringParams_.cosThbinEdge, steeringParams_.cosThbinWidth), LinearAxis(2011, 1));

  dataHist_.getAxis(0)->setLowerLimit(steeringParams_.minFitEnergy);
  dataHist_.getAxis(0)->setUpperLimit(steeringParams_.maxFitEnergy);
  dataHist_.getAxis(1)->setLowerLimit(steeringParams_.minCosth);
  dataHist_.getAxis(1)->setUpperLimit(steeringParams_.maxCosth);
  // fill in the histogram with the data
  bin(sample_, dataHist_, binner);
  data_histogram_constructed_=true;
}

void Sterilizer::ConstructSimulationHistogram(){
  if(not simulation_loaded_)
    throw std::runtime_error("No simulation has been loaded. Cannot construct simulation histogram.");
  if(not data_histogram_constructed_)
    throw std::runtime_error("Data histogram needs to be constructed before simulation histogram.");
  simHist_ = makeEmptyHistogramCopy(dataHist_);
  bin(mainSimulation_, simHist_, binner);
  simulation_histogram_constructed_=true;
}

/*************************************************************************************************************
 * Functions to obtain distributions
 * **********************************************************************************************************/

marray<double,3> Sterilizer::GetDataDistribution() const {
    if(not data_histogram_constructed_)
      throw std::runtime_error("Data histogram needs to be constructed before asking for it.");
    marray<double,3> array {static_cast<size_t>(dataHist_.getBinCount(2)),
                            static_cast<size_t>(dataHist_.getBinCount(1)),
                            static_cast<size_t>(dataHist_.getBinCount(0))};

    for(size_t iy=0; iy<dataHist_.getBinCount(2); iy++){ // year
      for(size_t ic=0; ic<dataHist_.getBinCount(1); ic++){ // zenith
        for(size_t ie=0; ie<dataHist_.getBinCount(0); ie++){ // energy
          auto itc = static_cast<likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(dataHist_(ie,ic,iy));
          array[iy][ic][ie] = itc.size();
        }
      }
    }
    return array;
}

marray<double,3> Sterilizer::GetExpectation(std::vector<double> nuisance) const {
  if(not simulation_histogram_constructed_)
    throw std::runtime_error("Simulation histogram needs to be constructed before asking for distributions.");

  marray<double,3> array {static_cast<size_t>(simHist_.getBinCount(2)),
                          static_cast<size_t>(simHist_.getBinCount(1)),
                          static_cast<size_t>(simHist_.getBinCount(0))};
  std::fill(array.begin(),array.end(),0);

  auto weighter = DFWM(nuisance);
  for(size_t iy=0; iy<simHist_.getBinCount(2); iy++){ // year
    for(size_t ic=0; ic<simHist_.getBinCount(1); ic++){ // zenith
      for(size_t ie=0; ie<simHist_.getBinCount(0); ie++){ // energy
        auto itc = static_cast<likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(simHist_(ie,ic,iy));
        double expectation=0;
        for(auto event : itc.entries()){
          expectation+=weighter(event);
        }
        array[iy][ic][ie] = expectation;
      }
    }
  }
  return array;
}

marray<double,3> Sterilizer::GetExpectation(Nuisance nuisance) const{
  return GetExpectation(ConvertNuisance(nuisance));
}

marray<double,3> Sterilizer::GetRealization(std::vector<double> nuisance, int seed) const{

  std::mt19937 rng;
  rng.seed(seed);

  std::cout << "construct weighter" << std::endl;
  auto weighter=DFWM(nuisance);

  double expected=0;
  std::vector<double> weights;
  for(const Event& e : mainSimulation_){
    auto w=weighter(e);
    if(std::isnan(w) || std::isinf(w) || w<0){
      std::cout << "Bad weight!" << std::endl;
      std::cout << e.cachedConvPionWeight  << ' ' << e.cachedConvKaonWeight << ' ' << e.cachedLivetime << ' ';
      std::cout << e.energy << ' ' << e.year << ' ' << w << std::endl;
    }
    weights.push_back(w);
    expected+=w;
  }

  std::vector<Event> realization= likelihood::generateSample(weights,mainSimulation_,expected,rng);
  auto realizationHist = makeEmptyHistogramCopy(dataHist_);
  bin(realization,realizationHist,binner);

  //if(realization.size() == 0 and not steeringParams_.quiet){
  if(realization.size() == 0){
    //std::cout << "No events generated. Expected events are "+std::to_string(expected) + "." << std::endl;
    throw std::runtime_error("No events generated. Expected events are "+std::to_string(expected));
  }

  marray<double,3> array {static_cast<size_t>(realizationHist.getBinCount(2)),
                          static_cast<size_t>(realizationHist.getBinCount(1)),
                          static_cast<size_t>(realizationHist.getBinCount(0))};
  std::fill(array.begin(),array.end(),0);

  for(size_t iy=0; iy<realizationHist.getBinCount(2); iy++){ // year
    for(size_t ic=0; ic<realizationHist.getBinCount(1); ic++){ // zenith
      for(size_t ie=0; ie<realizationHist.getBinCount(0); ie++){ // energy
        auto itc = static_cast<likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(realizationHist(ie,ic,iy));
        array[iy][ic][ie] = itc.size();
      }
    }
  }

  return array;
}

marray<double,3> Sterilizer::GetRealization(Nuisance nuisance, int seed) const {
  return GetRealization(ConvertNuisance(nuisance),seed);
}

/*************************************************************************************************************
 * Functions to construct likelihood problem and evaluate it
 * **********************************************************************************************************/

void Sterilizer::ConstructLikelihoodProblem(Priors pr, Nuisance nuisanceSeed, NuisanceFlag fixedParams){
  if(not data_histogram_constructed_)
    throw std::runtime_error("Data histogram needs to be constructed before likelihood problem can be formulated.");
  if(not simulation_histogram_constructed_)
    throw std::runtime_error("Simulation histogram needs to be constructed before likelihood problem can be formulated.");

  fixedParams_=fixedParams;
  fitSeed_=nuisanceSeed;

  UniformPrior  positivePrior(0.0,std::numeric_limits<double>::infinity());
  GaussianPrior normalizationPrior(pr.normCenter,pr.normWidth);
  GaussianPrior crSlopePrior(pr.crSlopeCenter,pr.crSlopeWidth);
  LimitedGaussianPrior simple_domEffPrior(pr.domEffCenter,pr.domEffWidth,-0.1,0.3);
  GaussianPrior kaonPrior(pr.piKRatioCenter,pr.piKRatioWidth);
  GaussianPrior nanPrior(pr.nuNubarRatioCenter,pr.nuNubarRatioWidth);
  GaussianPrior ZCPrior(0.0,pr.zenithCorrectionMultiplier*GetZenithCorrectionScale());

  auto llhpriors=makePriorSet(normalizationPrior,
			      positivePrior,
			      positivePrior,
			      crSlopePrior,
			      simple_domEffPrior,
			      kaonPrior,
			      nanPrior,
			      ZCPrior);

  auto fitseedvec = ConvertNuisance(nuisanceSeed);

  prob_ = std::make_shared<LType>(likelihood::makeLikelihoodProblem<std::reference_wrapper<const Event>, 3, 8>(
      dataHist_, {simHist_}, llhpriors, {1.0}, likelihood::simpleDataWeighter(), DFWM,
      likelihood::poissonLikelihood(), fitseedvec ));
  prob_->setEvaluationThreadCount(steeringParams_.evalThreads);

  likelihood_problem_constructed_=true;
}

// look up zenith correction scale for a particular flux model.
double Sterilizer::GetZenithCorrectionScale() const
{
  std::map<std::string,double> delta_alpha {
    {"HondaGaisser",8./7.},
    {"CombinedGHandHG_H3a_QGSJET",4. /7.},
    {"CombinedGHandHG_H3a_SIBYLL2",8./7.},
    {"PolyGonato_QGSJET-II-04",0.5},
    {"PolyGonato_SIBYLL2",1.0},
    {"ZatsepinSokolskaya_pamela_QGSJET",5./7.},
    {"ZatsepinSokolskaya_pamela_SIBYLL2",5./7.},
                  };
  if( delta_alpha.find(steeringParams_.modelName) == delta_alpha.end() )
    throw std::runtime_error("Jordi delta key not found. Aborting.");
  return delta_alpha[steeringParams_.modelName];
}

double Sterilizer::EvalLLH(std::vector<double> nuisance) const {
  if(not likelihood_problem_constructed_)
    throw std::runtime_error("Likelihood problem has not been constructed..");
  return -prob_->evaluateLikelihood(nuisance);
}

double Sterilizer::EvalLLH(Nuisance nuisance) const {
  return EvalLLH(ConvertNuisance(nuisance));
}

FitResult Sterilizer::MinLLH() const {
  if(not likelihood_problem_constructed_)
    throw std::runtime_error("Likelihood problem has not been constructed..");

  std::vector<double> seed=prob_->getSeed();
  std::vector<unsigned int> fixedIndices;

  std::vector<bool> FixVec=ConvertNuisanceFlag(fixedParams_);
  for(size_t i=0; i!=FixVec.size(); i++)
      if(FixVec[i]) fixedIndices.push_back(i);

  LBFGSB_Driver minimizer;

  minimizer.addParameter(seed[0],.001,0.0);
  minimizer.addParameter(seed[1],.001,0.0);
  minimizer.addParameter(seed[2],.01,0.0);
  minimizer.addParameter(seed[3],.005);
  minimizer.addParameter(seed[4],.005,-.1,.3);
  minimizer.addParameter(seed[5],.01,0.0);
  minimizer.addParameter(seed[6],.001,0.0,2.0);
  minimizer.addParameter(seed[7],.001,-1.0,1.0);

  minimizer.setChangeTolerance(1e-5);
  minimizer.setHistorySize(20);

  for(auto idx : fixedIndices)
    minimizer.fixParameter(idx);

  FitResult result;
  result.succeeded=DoFitLBFGSB(*prob_, minimizer);
  result.likelihood=minimizer.minimumValue();
  result.params=ConvertVecToNuisance(minimizer.minimumPosition());
  result.nEval=minimizer.numberOfEvaluations();
  result.nGrad=minimizer.numberOfEvaluations();

  return result;
}

/*************************************************************************************************************
 * Functions to change the sterile neutrino hypothesis
 * **********************************************************************************************************/

void Sterilizer::SetSterileNuParams(SterileNuParams snp){
  if(not simulation_loaded_)
    throw std::runtime_error("No simulation has been loaded. Cannot weight to sterile hypothesis without simulation.");
  if(not mc_generation_weighter_constructed_)
    throw std::runtime_error("MonteCarlo generation weighter has to be constructed first.");
  if(not xs_weighter_constructed_)
    throw std::runtime_error("Cross section weighter has to be constructed first.");

  sterileNuParams_=snp;
  ConstructFluxWeighter();
  ConstructLeptonWeighter();
  WeightMC();
}

/*************************************************************************************************************
 * Functions to set options in the class
 * **********************************************************************************************************/

// Check that the directories where files are mean to be exist
bool Sterilizer::CheckDataPaths(DataPaths dp) const
{
 CheckDataPath(dp.compact_file_path);
 CheckDataPath(dp.squids_files_path);
 CheckDataPath(dp.prompt_squids_files_path);
 CheckDataPath(dp.xs_spline_path);
 CheckDataPath(dp.data_path);
 CheckDataPath(dp.mc_path);
 CheckDataPath(dp.oversize_path);
 CheckDataPath(dp.domeff_spline_path);
 CheckDataPath(dp.flux_splines_path);
 return true;
}

 // Check a directory exists and throw a relevant error otherwise. 
 bool Sterilizer::CheckDataPath(std::string p) const
 {
   struct stat info;
   bool status=true;
   if(p!="")
     {
       if( stat(p.c_str(), &info) != 0 )
         {
           status=false;
           throw std::runtime_error("cannot access "+ p);
         }
       else if( !(info.st_mode & S_IFDIR) )
         {
           status=false;
           throw std::runtime_error("is not a directory: " +p);
         }
     }
   else{
     std::cout<<"Warning, there are unset paths in DataPaths. Check you want this."<<std::endl;
     return false;
   }
   return status;
 }

// Given a human readable nuisance parameter set, make a nuisance vector
std::vector<double> Sterilizer::ConvertNuisance(Nuisance ns) const {
  std::vector<double> nuis;
  nuis.push_back(ns.normalization);
  nuis.push_back(ns.astroFlux);
  nuis.push_back(ns.promptFlux);
  nuis.push_back(ns.crSlope);
  nuis.push_back(ns.domEfficiency);
  nuis.push_back(ns.piKRatio);
  nuis.push_back(ns.nuNubarRatio);
  nuis.push_back(ns.zenithCorrection);
  return nuis;
}

// Given a human readable flag set, make a bool vector
std::vector<bool> Sterilizer::ConvertNuisanceFlag(NuisanceFlag ns) const {
  std::vector<bool> nuis;
  nuis.push_back(ns.normalization);
  nuis.push_back(ns.astroFlux);
  nuis.push_back(ns.promptFlux);
  nuis.push_back(ns.crSlope);
  nuis.push_back(ns.domEfficiency);
  nuis.push_back(ns.piKRatio);
  nuis.push_back(ns.nuNubarRatio);
  nuis.push_back(ns.zenithCorrection);
  return nuis;
}

// And go back to human readable
Nuisance Sterilizer::ConvertVecToNuisance(std::vector<double> vecns) const {
  Nuisance ns;
  ns.normalization = vecns[0];
  ns.astroFlux = vecns[1];
  ns.promptFlux = vecns[2];
  ns.crSlope = vecns[3];
  ns.domEfficiency = vecns[4];
  ns.piKRatio = vecns[5];
  ns.nuNubarRatio = vecns[6];
  ns.zenithCorrection = vecns[7];
  return ns;
}

// Report back the status of object construction
void Sterilizer::ReportStatus() const {
  std::cout<< "Data loaded:                   " << CheckDataLoaded() <<std::endl;
  std::cout<< "Sim loaded:                    " << CheckSimulationLoaded()  <<std::endl;
  std::cout<< "Dom eff spline constructed:    " << CheckDOMEfficiencySplinesConstructed()<<std::endl;
  std::cout<< "XS weighter constructed:       " << CheckCrossSectionWeighterConstructed() <<std::endl;
  std::cout<< "Flux weighter constructed:     " << CheckFluxWeighterConstructed()<<std::endl;
  std::cout<< "Lepton weighter constructed:   " << CheckLeptonWeighterConstructed()<<std::endl;
  std::cout<< "Data histogram constructed:    " << CheckDataHistogramConstructed()<<std::endl;
  std::cout<< "Oversize weighter constructed: " << CheckOversizeWeighterConstructed()<<std::endl;
  std::cout<< "Sim histogram constructed:     " << CheckSimulationHistogramConstructed()<<std::endl;
  std::cout<< "LLH problem constructed:       " << CheckLikelihoodProblemConstruction()<<std::endl;
}

std::string Sterilizer::CheckedFilePath(std::string FilePath) const {
  if(!steeringParams_.quiet) std::cout<<"Reading a file from path "<<FilePath<<std::endl;
  try{
    std::ifstream thefile(FilePath);
    if(thefile.good())
      return FilePath;
    else
      throw std::runtime_error("File " + FilePath + " does not exist!");
  }
  catch(std::exception &re)
    {
      throw std::runtime_error("File " + FilePath + " does not exist!");
    }
}

/*************************************************************************************************************
 * Functions to spit out event distributions and swallow
 * **********************************************************************************************************/

double Sterilizer::Swallow(marray<double,2> Data)
{
  double TotalWeight=0;
  sample_.clear();
  for(size_t i=0; i!=Data.extent(0); ++i)
    {
      Event e;
      e.energy       = Data[i][0];
      e.zenith       = Data[i][1];
      e.year         = Data[i][2];
      e.cachedWeight = Data[i][3];
      TotalWeight+=Data[i][3];
      sample_.push_back(e);
    }
  // remaking data histogram
  if(!steeringParams_.quiet) std::cout<<"Remaking data hist" <<std::endl;
  ConstructDataHistogram();
  if(!steeringParams_.quiet) std::cout<<"ReMaking sim hist" <<std::endl;
  ConstructSimulationHistogram();
  // TO DO Improve this
  if(!steeringParams_.quiet) std::cout<<"Reconstrucing likelihood problem" <<std::endl;
  ConstructLikelihoodProblem(Priors(), Nuisance(),NuisanceFlag());
  return TotalWeight;
}

marray<double,2> Sterilizer::SpitData() const
{
  marray<double,2> ReturnVec { sample_.size(), 4} ;
  for(size_t i=0; i!=sample_.size(); ++i)
    {
      ReturnVec[i][0]=sample_[i].energy;
      ReturnVec[i][1]=sample_[i].zenith;
      ReturnVec[i][2]=sample_[i].year;
      ReturnVec[i][3]=sample_[i].cachedWeight;
    }
  return ReturnVec;
}


double Sterilizer::SetupDataChallenge(int seed, Nuisance nuisance)
{
  return Swallow(SpitRealization(ConvertNuisance(nuisance),seed));
}


marray<double,2> Sterilizer::SpitRealization( Nuisance nuisance, int seed) const
{
  return SpitRealization(ConvertNuisance(nuisance), seed);
}


marray<double,2> Sterilizer::SpitRealization( std::vector<double> nuisance, int seed) const
{
  std::mt19937 rng;
  rng.seed(seed);

  auto weighter=DFWM(nuisance);
  double expected=0;
  std::vector<double> weights;
  for(const Event& e : mainSimulation_){
    auto w=weighter(e);
    if(std::isnan(w) || std::isinf(w) || w<0){
      std::cout << "Bad weight!" << std::endl;
      std::cout << e.cachedConvPionWeight  << ' ' << e.cachedConvKaonWeight << ' ' << e.cachedLivetime << ' ';
      std::cout << e.energy << ' ' << e.year << ' ' << w << std::endl;
    }
    weights.push_back(w);
    expected+=w;
  }

  std::vector<Event> realization= likelihood::generateSample(weights,mainSimulation_,expected,rng); 

  marray<double,2> ReturnVec { realization.size(), 4} ;

  for(size_t i=0; i!=realization.size(); ++i)
    {
      ReturnVec[i][0]=realization[i].energy;
      ReturnVec[i][1]=realization[i].zenith;
      ReturnVec[i][2]=realization[i].year;
      ReturnVec[i][3]=1.;
    }
  return ReturnVec;

}

marray<double,2> Sterilizer::SpitExpectation( Nuisance nuisance) const
{
  return SpitExpectation(ConvertNuisance(nuisance));
}

marray<double,2> Sterilizer::SpitExpectation( std::vector<double> nuisance) const
{
  marray<double,2> ReturnVec { simHist_.getBinCount(2)*simHist_.getBinCount(1)*simHist_.getBinCount(0), 4} ;
  auto weighter = DFWM(nuisance);
  auto EnergyBins=GetEnergyBinsMC();
  auto ZenithBins=GetZenithBinsMC();
  unsigned int count=0;
  for(size_t iy=0; iy<simHist_.getBinCount(2); iy++){ // year
    for(size_t ic=0; ic<simHist_.getBinCount(1); ic++){ // zenith
      for(size_t ie=0; ie<simHist_.getBinCount(0); ie++){ // energy
        auto itc = static_cast<likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(simHist_(ie,ic,iy));
        double expectation=0;
        for(auto event : itc.entries()){
          expectation+=weighter(event);
        }
        ReturnVec[count][0]=(EnergyBins[ie]+EnergyBins[ie+1])/2.;
        ReturnVec[count][1]=(ZenithBins[ic]+ZenithBins[ic+1])/2.;
        ReturnVec[count][2]=steeringParams_.years[iy];
        ReturnVec[count][3]=expectation;
        ++count;
      }
    }
  }
  return ReturnVec;
}

bool Sterilizer::SetupAsimov(Nuisance nuisance)
{
  SetupAsimov(ConvertNuisance(nuisance));
  return true;
}

void Sterilizer::SetupAsimov(std::vector<double> nuisance)
{
  Swallow(SpitExpectation(nuisance));
}

/*************************************************************************************************************
 * Functions to get bin edges
 * **********************************************************************************************************/

// Given a histogram reference h, get bin edges in dimension dim
 std::vector<double> Sterilizer::PullBinEdges(int dim, const HistType& h) const{
   std::vector<double> edges_i(h.getBinCount(dim));
   for(unsigned int j=0; j<h.getBinCount(dim); j++)
     edges_i[j]=h.getBinEdge(dim,j);
   edges_i.push_back(h.getBinEdge(dim,h.getBinCount(dim)-1)+h.getBinWidth(dim,h.getBinCount(dim)-1));
   return edges_i;
 }

 std::vector<double> Sterilizer::GetEnergyBinsData() const{
   return PullBinEdges(0,dataHist_);
 }

 std::vector<double> Sterilizer::GetZenithBinsData() const{
   return PullBinEdges(1,dataHist_);
 }

 std::vector<double> Sterilizer::GetEnergyBinsMC() const{
   return PullBinEdges(0,simHist_);
 }

 std::vector<double> Sterilizer::GetZenithBinsMC() const{
   return PullBinEdges(1,simHist_);
 }


} // close namespace SterileSearch
