#include "SterileSearch.h"

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
    sample_=loadExperimentalData(dataPaths_.data_path,steeringParams_.useBurnSample);
  } catch(std::exception& ex){
    std::cerr << "Problem loading experimental data: " << ex.what() << std::endl;
  }
  if(!steeringParams_.quiet)
    std::cout << "Loaded " << sample_.size() << " experimental events" << std::endl;
  data_loaded_=true;
}

void Sterilizer::LoadMC(){
    bool loadTargeted=true;

    std::map<unsigned int,double> livetime;
    if(!steeringParams_.useBurnSample)
      livetime=steeringParams_.fullLivetime;
    else
      livetime=steeringParams_.burnSampleLivetime;

    std::vector<std::string> simSetsToLoad;
    simSetsToLoad.push_back(steeringParams_.simToLoad);
    try{
      loadSimulatedData(mainSimulation_,dataPaths_.mc_path,livetime,simInfo,simSetsToLoad,loadTargeted);
    } catch(std::exception& ex){
      std::cerr << "Problem loading simulated data: " << ex.what() << std::endl;
    }
    if(!steeringParams_.quiet)
      std::cout << "Loaded " << mainSimulation_.size() << " events in main simulation set" << std::endl;
    simulation_loaded_=true;
}

void Sterilizer::LoadCompact(){
  try{
    std::string original_data_file = dataPaths_.data_path+simInfo[steeringParams_.simToLoad].filename;
    std::string original_simulation_file = dataPaths_.data_path+simInfo[steeringParams_.simToLoad].filename;
    unsplatData(dataPaths_.compact_file_path+"/"+steeringParams_.simToLoad+"_compact_data.dat",
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
}

void Sterilizer::WriteCompact() const {
  try{
    std::string original_data_file = dataPaths_.data_path+simInfo[steeringParams_.simToLoad].filename;
    std::string original_simulation_file = dataPaths_.data_path+simInfo[steeringParams_.simToLoad].filename;
    splatData(dataPaths_.compact_file_path+"/"+steeringParams_.simToLoad+"_compact_data.dat",
	      getFileChecksum(original_data_file)+getFileChecksum(original_simulation_file),sample_,mainSimulation_);
  } catch(std::runtime_error& re){
    std::cerr << re.what() << std::endl;
    std::cerr << "Failed to save compact data" << std::endl;
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
  std::vector<unsigned int> years=steeringParams_.years;
  for(size_t year_index=0; year_index<steeringParams_.years.size(); year_index++){
    domEffConv_[year_index] = std::unique_ptr<Splinetable>(new Splinetable(dataPaths_.domeff_spline_path+"/conv_IC"+std::to_string(years[year_index])+".fits"));
    domEffPrompt_[year_index] = std::unique_ptr<Splinetable>(new Splinetable(dataPaths_.domeff_spline_path+"/prompt_IC"+std::to_string(years[year_index])+".fits"));
  }
  dom_efficiency_splines_constructed_=true;
}


/*************************************************************************************************************
 * Functions to construct weighters
 * **********************************************************************************************************/

void Sterilizer::ConstructCrossSectionWeighter(){
  xsw_ = std::make_shared<LW::CrossSectionFromSpline>(static_cast<std::string>(dataPaths_.xs_spline_path),steeringParams_.xs_model_name);
  xs_weighter_constructed_=true;
}

void Sterilizer::ConstructFluxWeighter(){
  std::string sterile_neutrino_model_identifier = GetSterileNeutrinoModelIdentifier(sterileNuParams_);

  if(steeringParams_.useFactorization){
    std::string oscillation_model = "noint_atmospheric_"+sterile_neutrino_model_identifier+".hdf5";
    std::string atmospheric_model_pion = steeringParams_.modelName + "_pion";
    std::string atmospheric_model_kaon = steeringParams_.modelName + "_kaon";

    fluxPion_ = std::make_shared<LW::FactorizedSQUIDSFlux>(dataPaths_.squids_files_path + oscillation_model,
                                                           dataPaths_.flux_splines_path+atmospheric_model_pion+"_neutrino_spline.fits",
                                                           dataPaths_.flux_splines_path+atmospheric_model_pion+"_antineutrino_spline.fits");
    fluxKaon_ = std::make_shared<LW::FactorizedSQUIDSFlux>(dataPaths_.squids_files_path+oscillation_model,
                                                           dataPaths_.flux_splines_path+atmospheric_model_kaon+"_neutrino_spline.fits",
                                                           dataPaths_.flux_splines_path+atmospheric_model_kaon+"_antineutrino_spline.fits");
  } else {
      std::string flux_pion_filename = "pion_atmospheric_"+sterile_neutrino_model_identifier;
      std::string flux_kaon_filename = "kaon_atmospheric_"+sterile_neutrino_model_identifier;
      flux_pion_filename+="_"+steeringParams_.modelName;
      flux_kaon_filename+="_"+steeringParams_.modelName;
      fluxKaon_ = std::make_shared<LW::SQUIDSFlux>(dataPaths_.squids_files_path + flux_kaon_filename + ".hdf5");
      fluxPion_ = std::make_shared<LW::SQUIDSFlux>(dataPaths_.squids_files_path + flux_pion_filename + ".hdf5");
  }

  fluxPrompt_ = std::make_shared<LW::SQUIDSFlux>(dataPaths_.prompt_squids_files_path + "prompt_atmospheric_0.000000_0.000000.hdf5");
  flux_weighter_constructed_=true;
}


void Sterilizer::ConstructMonteCarloGenerationWeighter(){
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

  PionFluxWeighter_ = LW::LeptonWeighter(fluxPion_,xsw_,mcw_);
  KaonFluxWeighter_ = LW::LeptonWeighter(fluxKaon_,xsw_,mcw_);
  PromptFluxWeighter_ = LW::LeptonWeighter(fluxPrompt_,xsw_,mcw_);
  lepton_weighter_constructed_=true;
}

void Sterilizer::ConstructOversizeWeighter(){
  osw_=OversizeWeighter(dataPaths_.oversize_path+"/"+steeringParams_.oversizeFunction+".dat");
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
  initializeSimulationWeights(mainSimulation_,PionFluxWeighter_,KaonFluxWeighter_,PromptFluxWeighter_,osw_);
  simulation_initialized_=true;
}

/*************************************************************************************************************
 * Functions to construct histograms
 * **********************************************************************************************************/

void Sterilizer::ConstructDataHistogram(){
  if(not data_loaded_)
    throw std::runtime_error("No data has been loaded. Cannot construct data histogram.");

  dataHist_ = HistType(LogarithmicAxis(0, 0.1), LinearAxis(0, 0.1), LinearAxis(2010, 1));

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

  std::deque<Event> realization= likelihood::generateSample(weights,mainSimulation_,expected,rng);
  auto realizationHist = makeEmptyHistogramCopy(dataHist_);
  bin(realization,realizationHist,binner);

  marray<double,3> array {static_cast<size_t>(realizationHist.getBinCount(2)),
      static_cast<size_t>(realizationHist.getBinCount(1)),
      static_cast<size_t>(realizationHist.getBinCount(0))};

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

void Sterilizer::ConstructLikelihoodProblem(Priors priors, Nuisance nuisanceSeed){
  if(not data_histogram_constructed_)
    throw std::runtime_error("Data histogram needs to be constructed before likelihood problem can be formulated.");
  if(not simulation_histogram_constructed_)
    throw std::runtime_error("Simulation histogram needs to be constructed before likelihood problem can be formulated.");

  auto llhpriors = ConvertPriorSet(priors);
  auto fitseed   = ConvertNuisance(nuisanceSeed);

  prob_ = likelihood::makeLikelihoodProblem<std::reference_wrapper<const Event>, 3, 6>(
      dataHist_, {simHist_}, llhpriors, {1.0}, likelihood::simpleDataWeighter(), DFWM,
      likelihood::poissonLikelihood(), fitseed );
  prob_.setEvaluationThreadCount(steeringParams_.evalThreads);

  likelihood_problem_constructed_=true;
}

double Sterilizer::EvalLLH(std::vector<double> nuisance) const {
  if(not likelihood_problem_constructed_)
    throw std::runtime_error("Likelihood problem has not been constructed..");
  return -prob_.evaluateLikelihood(nuisance);
}

double Sterilizer::EvalLL(Nuisance nuisance) const {
  return EvalLLH(ConvertNuisance(nuisance));
}

FitResult Sterilizer::MinLLH(NuisanceFlag fixedParams) const {
  if(not likelihood_problem_constructed_)
    throw std::runtime_error("Likelihood problem has not been constructed..");

  std::vector<double> seed=prob_.getSeed();
  std::vector<unsigned int> fixedIndices;

  std::vector<bool> FixVec=ConvertNuisanceFlag(fixedParams);
  for(size_t i=0; i!=FixVec.size(); ++i)
      if(FixVec[i]) fixedIndices.push_back(i);

  return DoFitLBFGSB(prob_, seed, fixedIndices);
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
}

 // Check a directory exists and throw a relevant error otherwise. 
 bool Sterilizer::CheckDataPath(std::string p) const
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
     std::cout<<"Warning, there are unset paths in DataPaths. Check you want this."<<std::endl;
 }


// Given a human readable prior set, make a weaverized version
template<typename... PriorTypes>
    FixedSizePriorSet<PriorTypes...> Sterilizer::ConvertPriorSet(Priors pr) const
{
  // construct continuous nuisance priors      
  UniformPrior  positivePrior(0.0,std::numeric_limits<double>::infinity());
  GaussianPrior normalizationPrior(pr.normCenter,pr.normWidth);
  GaussianPrior crSlopePrior(pr.crSlopeCenter,pr.crSlopeWidth);
  UniformPrior  simple_domEffPrior(pr.domEffCenter,pr.domEffWidth);
  GaussianPrior kaonPrior(pr.piKRatioCenter,pr.piKRatioWidth);
  GaussianPrior nanPrior(pr.nuNubarRatioCenter,pr.nuNubarRatioWidth);

  // construct zenith correction prior                   
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
  double alpha = delta_alpha[steeringParams_.modelName];

  GaussianPrior ZCPrior(0.0,pr.zenithCorrectionMultiplier*alpha);

  // make and return priorset  
  return makePriorSet(normalizationPrior,
                             positivePrior, 
                             positivePrior,
                             crSlopePrior,
                             simple_domEffPrior,
                             kaonPrior,
                             nanPrior,
                             ZCPrior);
}

// Given a human readable nuisance parameter set, make a nuisance vector
std::vector<double> Sterilizer::ConvertNuisance(Nuisance ns) const {
  std::vector<double> nuis;
  nuis.append(ns.normalization);
  nuis.append(ns.astroFlux);
  nuis.append(ns.promptFlux);
  nuis.append(ns.crSlope);
  nuis.append(ns.domEfficiency);
  nuis.append(ns.piKRatio);
  nuis.append(ns.nuNubarRatio);
  nuis.append(ns.zenithCorrection);
  return nuis;
}

// Given a human readable flag set, make a bool vector
std::vector<bool> Sterilizer::ConvertNuisanceFlag(NuisanceFlag ns) const {
  return std::vector<bool> nuis{
    ns.normalization,
      ns.astroFlux,
      ns.promptFlux,
      ns.crSlope,
      ns.domEfficiency,
      ns.piKRatio,
      ns.nuNubarRatio,
      ns.zenithCorrection
      }
}

// And go back to human readable
std::vector<double> Sterilizer::ConvertVecToNuisance(std::vector<double> vecns) const {
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

// Do the fit business
template<typename LikelihoodType>
FitResult Sterilizer::DoFitLBFGSB(LikelihoodType& likelihood, const std::vector<double>& seed,
		      std::vector<unsigned int> indicesToFix={}){
  using namespace likelihood;

  LBFGSB_Driver minimizer;
  minimizer.addParameter(seed[0],.001,0.0);
  minimizer.addParameter(seed[1],.001,0.0);
  minimizer.addParameter(seed[2],.01,0.0);
  minimizer.addParameter(seed[3],.005);
  minimizer.addParameter(seed[4],.005,-.1,.3);
  minimizer.addParameter(seed[5],.01,0.0);
  minimizer.addParameter(seed[6],.001,0.0,2.0);
  minimizer.addParameter(seed[7],.001,-1.0,1.0);

  for(auto idx : indicesToFix)
    minimizer.fixParameter(idx);

  minimizer.setChangeTolerance(1e-5);
  minimizer.setHistorySize(20);
  FitResult result;
  result.succeeded=minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));
  result.likelihood=minimizer.minimumValue();
  result.params=ConvertVecToNuisance(minimizer.minimumPosition());
  result.nEval=minimizer.numberOfEvaluations();
  result.nGrad=minimizer.numberOfEvaluations(); //gradient is always eval'ed with function

  return(result);
}
