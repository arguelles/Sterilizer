#include "SterileSearch.h"
#include "GenerationSpecifications.h"

namespace SterileSearch {

/*************************************************************************************************************
 * Constructor
 * **********************************************************************************************************/

Sterilizer::Sterilizer(DataPaths dataPaths, SteeringParams steeringParams, SterileNuParams snp):
  steeringParams_(steeringParams),dataPaths_(dataPaths),sterileNuParams_(snp),
  nus_atm_pion_(nuSQUIDSAtm<>(linspace(-1.,0.2,40),logspace(1.e2*units.GeV,1.e6*units.GeV,150),numneu,both,true)),
  nus_atm_kaon_(nuSQUIDSAtm<>(linspace(-1.,0.2,40),logspace(1.e2*units.GeV,1.e6*units.GeV,150),numneu,both,true)){

  if(!steeringParams_.quiet) std::cout<<"Sterilizer constructor: checking paths" <<std::endl;
  CheckDataPaths(dataPaths_);

  if(!steeringParams_.quiet) std::cout<<"Loading Splines" <<std::endl;
  LoadDOMEfficiencySplines();

  if(steeringParams_.readCompact){
    if(!steeringParams_.quiet) std::cout<<"Loading compact data" <<std::endl;
    if(steeringParams_.fastMode)
      LoadFastCompact();
    else
      LoadCompact();
  } else {
    if(!steeringParams_.quiet) std::cout<<"Loading data" <<std::endl;
    LoadData();
    if(!steeringParams_.quiet) std::cout<<"Loading MC" <<std::endl;
    LoadMC();
  }
  if(!(steeringParams_.readCompact && steeringParams_.fastMode))
    {
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
    }
  if(!steeringParams_.quiet) std::cout<<"Making data hist" <<std::endl;
  ConstructDataHistogram();
  if(!steeringParams_.quiet) std::cout<<"Making sim hist" <<std::endl;
  ConstructSimulationHistogram();
  if(!steeringParams_.quiet) std::cout<<"Construcing likelihood problem with default settings" <<std::endl;
  if(steeringParams_.fastMode)
    if(!steeringParams_.readCompact)
      SetupFastMode();
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

    try{
      auto simAction=[&](RecordID id, Event& e, int simYear, const simpleEffRate<Event>& domEff){
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
	domEffObject_=new simpleEffRate<Event>(domEffConv_[simYear].get(),setInfo.unshadowedFraction);
	auto callback=[&,simYear](RecordID id, Event& e){ simAction(id,e,simYear,*domEffObject_); };
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

void Sterilizer::LoadMCFromTextFile(){
  if (!steeringParams_.quiet)
    std::cout << "Loading MC from TXT file data." << std::endl;

  auto mc_dump= nusquids::quickread(CheckedFilePath(dataPaths_.data_path+"mc_dump.dat"));
  mainSimulation_.clear();

  for (unsigned int irow = 0; irow < mc_dump.extent(0); irow++){
    Event evt;

    if(mc_dump[irow][0] == 14)
      evt.primaryType = particleType::NuMu;
    else if(mc_dump[irow][0] == -14)
      evt.primaryType = particleType::NuMuBar;
    else
      throw std::runtime_error( "what particle" );
    evt.injectedEnergy=mc_dump[irow][1]; // neutrino energy
    evt.injectedMuonZenith=mc_dump[irow][2]; // neutrino direction
    evt.energy=mc_dump[irow][3]; // reconstructed energy
    evt.zenith=mc_dump[irow][4]; // reconstructed zenith
    evt.cachedWeight=mc_dump[irow][5]; // one event is one event
    evt.year=static_cast<unsigned int>(2011); // year of the event

    mainSimulation_.push_back(evt); // push it
  }

  simulation_loaded_=true;
  simplified_simulation_=true;
}

void Sterilizer::LoadCompact(){
  std::map<std::string,run> simInfo=GetSimInfo(dataPaths_.mc_path);
  try{
    auto simulation_information=simInfo.find(steeringParams_.simToLoad);
    if(simulation_information==simInfo.end())
      throw std::runtime_error("Could not find " + steeringParams_.simToLoad + " in simulation list");

    std::string original_data_file = CheckedFilePath(dataPaths_.data_path+"IC86.h5");
    std::string original_simulation_file = CheckedFilePath(dataPaths_.mc_path+(*simulation_information).second.filename);

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


void Sterilizer::LoadFastCompact(){
  try{
    unsplatData(CheckedFilePath(dataPaths_.compact_file_path+"/fast_sim_hist_"+std::to_string(sterileNuParams_.modelId)+".dat"),
        0,sample_,mainSimulation_);

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



bool Sterilizer::WriteFastCompact() const {
  std::map<std::string,run> simInfo=GetSimInfo(dataPaths_.mc_path);
  try{
    splatData(dataPaths_.compact_file_path+"/fast_sim_hist_"+std::to_string(sterileNuParams_.modelId)+".dat",
	      0,sample_,auxSimulation_);
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
  }
  DFWM.SetSplines(domEffConv_);
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

    if(dataPaths_.use_simple_filename){
      flux_pion_filename = "pion_"+std::to_string(sterileNuParams_.modelId);
      flux_kaon_filename = "kaon_"+std::to_string(sterileNuParams_.modelId);
    }

    if(steeringParams_.calculate_nusquids_on_the_fly){
      ConstructNuSQuIDSObjects();
      fluxKaon_ = std::make_shared<LW::SQUIDSFlux>(std::move(nus_atm_kaon_));
      fluxPion_ = std::make_shared<LW::SQUIDSFlux>(std::move(nus_atm_pion_));
    } else {
      fluxKaon_ = std::make_shared<LW::SQUIDSFlux>(CheckedFilePath(dataPaths_.squids_files_path + flux_kaon_filename + ".hdf5"));
      fluxPion_ = std::make_shared<LW::SQUIDSFlux>(CheckedFilePath(dataPaths_.squids_files_path + flux_pion_filename + ".hdf5"));
    }
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
      //     std::cout<< e.injectedEnergy<<", " <<e.injectedMuonZenith<<", " <<e.cachedConvPionWeight<<", " <<e.cachedConvKaonWeight<<std::endl;
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

  if(steeringParams_.fastMode && steeringParams_.readCompact){
    std::cout<<"This functionality is not available in fast / compact mode!"<<std::endl;
    assert(0);
  }

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
			      positivePrior,
			      kaonPrior,
			      nanPrior,
			      ZCPrior);

  auto fitseedvec = ConvertNuisance(nuisanceSeed);




  prob_ = std::make_shared<LType>(likelihood::makeLikelihoodProblem<std::reference_wrapper<const Event>, 3, 8>(
													       dataHist_, {simHist_}, llhpriors, {1.0}, simpleLocalDataWeighter(), DFWM,
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

  minimizer.addParameter(seed[0],.001,0.5,3);
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


double Sterilizer::SetupDataChallenge(int seed, Nuisance nuisance_dc, SterileNuParams snp_dc)
{
  // Save these two things for later, since we'll be adjusting them
  bool my_quiet=steeringParams_.quiet;
  SterileNuParams my_snp=sterileNuParams_;

  // Set the parameter point to the data challenge point and make realization
  steeringParams_.quiet=true;
  SetSterileNuParams(snp_dc);
  marray<double,2> TheRealization=SpitRealization(ConvertNuisance(nuisance_dc),seed);

  // Set the parameter point back to the original one
  SetSterileNuParams(my_snp);
  steeringParams_.quiet=my_quiet;

  // Set the realization as data
  return Swallow(TheRealization);
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


// Warning - do not call this publically (runs in constructor)
//  Setup fast mode (one sim event per bin)
void Sterilizer::SetupFastMode()
{
  std::cout<<"Setting up fast mode"<<std::endl;
  marray<double,2> ReturnVec { simHist_.getBinCount(2)*simHist_.getBinCount(1)*simHist_.getBinCount(0), 4} ;
  auto weighter = DFWM(ConvertNuisance(Nuisance()));
  auto EnergyBins=GetEnergyBinsMC();
  auto ZenithBins=GetZenithBinsMC();
  unsigned int count=0;
  auxSimulation_.clear();

  for(size_t iy=0; iy<simHist_.getBinCount(2); iy++){ // year
    for(size_t ic=0; ic<simHist_.getBinCount(1); ic++){ // zenith
      for(size_t ie=0; ie<simHist_.getBinCount(0); ie++){ // energy
        auto itc = static_cast<likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(simHist_(ie,ic,iy));
        double expectationnu=0;
        double expectationnubar=0;
	Event evnu, evnubar;

	evnu.leptonEnergyFraction=0;
	evnu.totalColumnDepth=0;
	evnu.inelasticityProbability=0;
	evnu.intX=0;
	evnu.intY=0;
	evnu.cutL3=false;
	evnu.paraboloidStatus=0;
	evnu.cachedConvPionWeight=0;
	evnu.cachedConvKaonWeight=0;
	evnu.cachedPromptWeight=0;
	evnu.cachedWeight=0;
	evnu.injectedEnergy =0;
	evnu.energy =0;
	evnu.zenith =0;
	evnu.livetime =0;
	evnu.primaryType=particleType::MuMinus;

	evnubar.leptonEnergyFraction=0;
	evnubar.totalColumnDepth=0;
	evnubar.inelasticityProbability=0;
	evnubar.intX=0;
	evnubar.intY=0;
	evnubar.cutL3=false;
	evnubar.paraboloidStatus=0;	
	evnubar.cachedConvPionWeight=0;
	evnubar.cachedConvKaonWeight=0;
	evnubar.cachedPromptWeight=0;
	evnubar.cachedWeight=0;
	evnubar.injectedEnergy =0;
	evnubar.energy =0;
	evnubar.zenith =0;
	evnubar.livetime =0;
	evnubar.primaryType=particleType::MuPlus;

        for(auto event : itc.entries()){
          double weight=weighter(event);
	  //	  std::cout<<"This evt weight:" << weight<<std::endl;
	  Event theev(event);
	  
	  if(theev.primaryType==particleType::MuMinus)
	    {
	      //	      std::cout<<"nu"<<std::endl;
	      expectationnu+=weight;
	
	      evnu.injectedEnergy += theev.injectedEnergy*weight;
	      evnu.energy   += theev.energy*weight;
	      evnu.zenith   += theev.zenith*weight;
	      evnu.livetime += theev.livetime*weight;
	      
	      evnu.cachedConvPionWeight += theev.cachedConvPionWeight;
	      evnu.cachedConvKaonWeight += theev.cachedConvKaonWeight;
	      evnu.cachedPromptWeight   += theev.cachedPromptWeight;
	      evnu.cachedWeight         += theev.cachedWeight;
	      evnu.year=theev.year;
	 	
	    }
	  else if (theev.primaryType==particleType::MuPlus)
	    {
	      //     std::cout<<"nubar"<<std::endl;
	      expectationnubar+=weight;
	
	      evnubar.injectedEnergy += theev.injectedEnergy*weight;
	      evnubar.energy   += theev.energy*weight;
	      evnubar.zenith   += theev.zenith*weight;
	      evnubar.livetime += theev.livetime*weight;
	      
	      evnubar.cachedConvPionWeight += theev.cachedConvPionWeight;
	      evnubar.cachedConvKaonWeight += theev.cachedConvKaonWeight;
	      evnubar.cachedPromptWeight   += theev.cachedPromptWeight;
	      evnubar.cachedWeight         += theev.cachedWeight;
	      evnubar.year=theev.year;

	    }
	}
	if(expectationnu>0)
	  {
	    evnu.energy/=expectationnu;
	    evnu.zenith/=expectationnu;
	    evnu.injectedEnergy/=expectationnu;
	    evnu.livetime/=expectationnu;
	    domEffObject_->setCache(evnu);
	    auxSimulation_.push_back(evnu);
	  }
	if(expectationnubar>0)
	  {
	    evnubar.energy/=expectationnubar;
	    evnubar.zenith/=expectationnubar;
	    evnubar.injectedEnergy/=expectationnubar;
	    evnubar.livetime/=expectationnubar;
	    domEffObject_->setCache(evnubar);
	    auxSimulation_.push_back(evnubar);    
	  }
		
      }
    }
  }
  simHist_=makeEmptyHistogramCopy(dataHist_);
  bin(auxSimulation_, simHist_, binner);
  fastmode_constructed_=true;
  std::cout<<"done constructing fast mode"<<std::endl;
  
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

double Sterilizer::SetupAsimovForAlternativeHypothesis(Nuisance nuisance, SterileNuParams snp_dc) {
  return SetupAsimovForAlternativeHypothesis(ConvertNuisance(nuisance),snp_dc);
}

double Sterilizer::SetupAsimovForAlternativeHypothesis(std::vector<double> nuisance_dc, SterileNuParams snp_dc) {
  // Save these two things for later, since we'll be adjusting them
  bool my_quiet=steeringParams_.quiet;
  SterileNuParams my_snp=sterileNuParams_;

  // Set the parameter point to the data challenge point and make realization
  steeringParams_.quiet=true;
  SetSterileNuParams(snp_dc);
  marray<double,2> TheAsimovRealization=SpitExpectation(nuisance_dc);

  // Set the parameter point back to the original one
  SetSterileNuParams(my_snp);
  steeringParams_.quiet=my_quiet;

  // Set the realization as data
  return Swallow(TheAsimovRealization);
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

/*************************************************************************************************************
 * Functions to construct nusquids state on the fly
 * **********************************************************************************************************/

void Sterilizer::ConstructNuSQuIDSObjects(){
  using namespace nusquids;

  nus_atm_kaon_.Set_TauRegeneration(false);
  nus_atm_pion_.Set_TauRegeneration(false);

  nus_atm_kaon_.Set_MixingAngle(0,1,0.563942);
  nus_atm_kaon_.Set_MixingAngle(0,2,0.154085);
  nus_atm_kaon_.Set_MixingAngle(1,2,0.785398);
  nus_atm_kaon_.Set_MixingAngle(0,3,sterileNuParams_.th14);
  nus_atm_kaon_.Set_MixingAngle(1,3,sterileNuParams_.th24);
  nus_atm_kaon_.Set_MixingAngle(2,3,sterileNuParams_.th34);

  nus_atm_kaon_.Set_SquareMassDifference(1,7.65e-05);
  nus_atm_kaon_.Set_SquareMassDifference(2,0.00247);
  nus_atm_kaon_.Set_SquareMassDifference(3,sterileNuParams_.dm41sq);

  nus_atm_kaon_.Set_CPPhase(0,2,0.0);
  nus_atm_kaon_.Set_CPPhase(0,3,sterileNuParams_.del14);
  nus_atm_kaon_.Set_CPPhase(1,3,sterileNuParams_.del24);

  nus_atm_pion_.Set_MixingAngle(0,1,0.563942);
  nus_atm_pion_.Set_MixingAngle(0,2,0.154085);
  nus_atm_pion_.Set_MixingAngle(1,2,0.785398);
  nus_atm_pion_.Set_MixingAngle(0,3,sterileNuParams_.th14);
  nus_atm_pion_.Set_MixingAngle(1,3,sterileNuParams_.th24);
  nus_atm_pion_.Set_MixingAngle(2,3,sterileNuParams_.th34);

  nus_atm_pion_.Set_SquareMassDifference(1,7.65e-05);
  nus_atm_pion_.Set_SquareMassDifference(2,0.00247);
  nus_atm_pion_.Set_SquareMassDifference(3,sterileNuParams_.dm41sq);
  nus_atm_pion_.Set_CPPhase(0,2,0.0);
  nus_atm_pion_.Set_CPPhase(0,3,sterileNuParams_.del14);
  nus_atm_pion_.Set_CPPhase(1,3,sterileNuParams_.del24);

  double error = 1.0e-10;
  // setup integration settings
  nus_atm_pion_.Set_GSL_step(gsl_odeiv2_step_rkf45);
  //nus_atm_pion.Set_GSL_step(gsl_odeiv2_step_rk4);
  nus_atm_pion_.Set_rel_error(error);
  nus_atm_pion_.Set_abs_error(error);

  nus_atm_kaon_.Set_GSL_step(gsl_odeiv2_step_rkf45);
  //nus_atm_kaon_.Set_GSL_step(gsl_odeiv2_step_rk4);
  nus_atm_kaon_.Set_rel_error(error);
  nus_atm_kaon_.Set_abs_error(error);

  if(not steeringParams_.quiet){
    nus_atm_kaon_.Set_ProgressBar(true);
    nus_atm_pion_.Set_ProgressBar(true);
  } else {
    nus_atm_kaon_.Set_ProgressBar(false);
    nus_atm_pion_.Set_ProgressBar(false);
  }

  // read file
  marray<double,2> input_pion_flux = quickread(CheckedFilePath(dataPaths_.initial_flux_files_path + "/" + "initial_pion_atmopheric_" + steeringParams_.modelName + ".dat"));
  marray<double,2> input_kaon_flux = quickread(CheckedFilePath(dataPaths_.initial_flux_files_path + "/" + "initial_kaon_atmopheric_" + steeringParams_.modelName + ".dat"));

  marray<double,4> inistate_kaon {nus_atm_kaon_.GetNumCos(),nus_atm_kaon_.GetNumE(),2,numneu};
  std::fill(inistate_kaon.begin(),inistate_kaon.end(),0);

  marray<double,1> cos_range = nus_atm_kaon_.GetCosthRange();
  marray<double,1> e_range = nus_atm_kaon_.GetERange();
  for ( int ci = 0 ; ci < nus_atm_kaon_.GetNumCos(); ci++){
    for ( int ei = 0 ; ei < nus_atm_kaon_.GetNumE(); ei++){
      double enu = e_range[ei]/units.GeV;
      double cth = cos_range[ci];

      inistate_kaon[ci][ei][0][0] = 0.;
      inistate_kaon[ci][ei][0][1] = input_kaon_flux[ci*e_range.size() + ei][2];
      inistate_kaon[ci][ei][0][2] = 0.;
      inistate_kaon[ci][ei][0][3] = 0.;

      inistate_kaon[ci][ei][1][0] = 0.;
      inistate_kaon[ci][ei][1][1] = input_kaon_flux[ci*e_range.size() + ei][3];
      inistate_kaon[ci][ei][1][2] = 0.;
      inistate_kaon[ci][ei][1][3] = 0.;
    }
  }

  nus_atm_kaon_.Set_initial_state(inistate_kaon,flavor);
  nus_atm_kaon_.EvolveState();

  marray<double,4> inistate_pion {nus_atm_pion_.GetNumCos(),nus_atm_pion_.GetNumE(),2,numneu};
  std::fill(inistate_pion.begin(),inistate_pion.end(),0);

  for ( int ci = 0 ; ci < nus_atm_pion_.GetNumCos(); ci++){
    for ( int ei = 0 ; ei < nus_atm_pion_.GetNumE(); ei++){
      double enu = e_range[ei]/units.GeV;
      double cth = cos_range[ci];

      inistate_pion[ci][ei][0][0] = 0.;
      inistate_pion[ci][ei][0][1] = input_pion_flux[ci*e_range.size() + ei][2];
      inistate_pion[ci][ei][0][2] = 0.;
      inistate_pion[ci][ei][0][3] = 0.;

      inistate_pion[ci][ei][1][0] = 0.;
      inistate_pion[ci][ei][1][1] = input_pion_flux[ci*e_range.size() + ei][3];
      inistate_pion[ci][ei][1][2] = 0.;
      inistate_pion[ci][ei][1][3] = 0.;
    }
  }

  nus_atm_pion_.Set_initial_state(inistate_pion,flavor);
  nus_atm_pion_.EvolveState();
}

} // close namespace SterileSearch
