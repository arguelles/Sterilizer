#include "SterileSearch.h"


/*************************************************************************************************************
 * Implementation auxiliary functions
 * **********************************************************************************************************/

std::string GetSterileNeutrinoModelIdentifier(SterileNeutrinoParameters snp){
  return snp.modelId+"_"+std::to_string(snp.dm41sq)+"_"+std::to_string(snp.th14)+"_"+std::to_string(snp.th24)+"_"+std::to_string(snp.th34)+"_"+std::to_string(snp.del14)+"_"+std::to_string(snp.del24);
}

/*************************************************************************************************************
 * Functions to read and write data
 * **********************************************************************************************************/

void Sterilizer::LoadData(){
  try{
    sample_=loadExperimentalData(dataPaths.data_path,steeringParams_.useBurnsample);
  } catch(std::exception& ex){
    std::cerr << "Problem loading experimental data: " << ex.what() << std::endl;
    return(1);
  }
  if(!quiet)
    std::cout << "Loaded " << sample_.size() << " experimental events" << std::endl;
  data_loaded=true;
}

void Sterilizer::LoadMC(){
    bool loadTargeted=true;
    std::vector<std::string> simSetsToLoad;
    simSetsToLoad.push_back(steeringParams_.simToLoad_);
    try{
      loadSimulatedData(mainSimulation_,dataPaths.mc_path,livetime,simInfo,simSetsToLoad,loadTargeted);
    } catch(std::exception& ex){
      std::cerr << "Problem loading simulated data: " << ex.what() << std::endl;
      return(1);
    }
    if(!quiet)
      std::cout << "Loaded " << mainSimulation_.size() << " events in main simulation set" << std::endl;
    simulation_loaded=true;
}

void Sterilizer::LoadCompact(){
  try{
    unsplatData(dataPaths_.compact_data_path+"/"+steeringParams_.simToLoad+"_compact_data.dat",getFileChecksum(argv[0]),sample_,mainSimulation_);
    if(!quiet){
      std::cout << "Loaded " << sample_.size() << " experimental events." << std::endl;
      std::cout << "Loaded " << mainSimulation_.size() << " events in main simulation set." << std::endl;
    }
  } catch(std::runtime_error& re){
    std::cerr << re.what() << std::endl;
    std::cerr << "Failed to load compact data" << std::endl;
    return(1);
  }
  data_loaded=true;
  simulation_loaded=true;
}

void Sterilizer::WriteCompact() const {
  try{
    splatData(dataPaths_.compact_data_path+"/"+steeringParams_.simToLoad+"_compact_data.dat",
	      getFileChecksum(argv[0]),sample_,mainSimulation_);
  } catch(std::runtime_error& re){
    std::cerr << re.what() << std::endl;
    std::cerr << "Failed to save compact data" << std::endl;
    return(1);
  }
}

void Sterilizer::ClearData(){
  sample_.clear();
}

void Sterilizer::ClearSimulation(){
  mainSimulation_.clear();
}

bool Sterilizer::CheckDataLoaded() const {
  return(data_loaded);
}

bool Sterilizer:;CheckSimulationLoaded() const {
  return(simulation_to_loaded);
}

/*************************************************************************************************************
 * Functions to load to load DOM efficiency splines
 * **********************************************************************************************************/

void Sterilizer::LoadDOMEfficiencySplines(std::vector<unsigned int> years){
  for(size_t year_index=0; year_index<years.size(); year_index++){
    domEffConv_[year_index] = std::unique_ptr<Splinetable>(new Splinetable(dataPaths_.domeff_spline_path+"/conv_IC"+std::to_string(years[year_index]+".fits"));
    domEffPrompt_[year_index] = std::unique_ptr<Splinetable>(new Splinetable(dataPaths_.domeff_spline_path+"/prompt_IC"+std::to_string(years[year_index])+".fits"));
  }
  dom_efficiency_splines_loaded_=true;
}

bool Sterilizer::CheckDOMEfficiencySplinesLoaded() const {
  return(dom_efficiency_splines_loaded_);
}

/*************************************************************************************************************
 * Functions to construct weighters
 * **********************************************************************************************************/

void Sterilizer::ConstructCrossSectionWeighter(){
  xsw_ = std::make_shared<LW::CrossSectionFromSpline>(static_cast<std::string>(dataPaths.xs_spline_path),steeringParams_.xs_model_name);
  cross_section_weighter_constructed_=true;
}

bool SterileSearch::CheckCrossSectionWeighterConstructed() const {
  return(cross_section_weighter_constructed_);
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
  } else{
      std::string flux_pion_filename = "pion_atmospheric_"+sterile_neutrino_model_identifier;
      std::string flux_kaon_filename = "kaon_atmospheric_"+sterile_neutrino_model_identifier;
      if(model_name != ""){
        flux_pion_filename+="_"+model_name;
        flux_kaon_filename+="_"+model_name;
      }
      fluxKaon_ = std::make_shared<LW::SQUIDSFlux>(dataPaths_.squids_files_path + flux_kaon_filename + ".hdf5");
      fluxPion_ = std::make_shared<LW::SQUIDSFlux>(dataPaths_.squids_files_path + flux_pion_filename + ".hdf5");
  }

  fluxPrompt_ = std::make_shared<LW::SQUIDSFlux>(dataPaths_.prompt_squids_files_path + "prompt_atmospheric_0.000000_0.000000.hdf5");
  flux_weighter_constructed_=true;
}

bool Sterilizer::CheckFluxWeighterConstructed() const {
  return(flux_weighter_constructed_);
}

void Sterilizer::ConstructMonteCarloGenerationWeighter(){
  std::vector<std::string> simSetsToLoad;
  simSetsToLoad.push_back(steeringParams_.simToLoad);
  for( std::string sim_name : simSetsToLoad )
    mcw_.addGenerationSpectrum(simInfo.find(sim_name)->second.details);
  mc_generation_weighter_constructed_=true;
}

bool Sterilizer::CheckFluxWeighterConstructed() const {
  return(flux_weighter_constructed_);
}

void Sterilizer::ConstructLeptonWeighter(){
  if(not mc_generation_weighter_constructed_)
    throw std::runtime_error("MonteCarlo generation weighter has to be constructed first.");
  if(not flux_weighter_constructed_)
    throw std::runtime_error("Flux weighter has to be constructed first.");
  if(not cross_section_weighter_constructed_)
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

bool SterileSearch::CheckLeptonWeighterConstructed() const {
  return(lepton_weighter_constructed_);
}

/*************************************************************************************************************
 * Functions to initialize the MC weights
 * **********************************************************************************************************/

void Sterilizer::WeightMC(){
  if(not lepton_weighter_constructed_)
    throw std::runtime_error("LeptonWeighter has to be constructed first.");
  if(not oversize_weighter_constructed_)
    throw std::runtime_error("OversizeWeighter has to be constructed first.");
  if(not simulation_loaded))
    throw std::runtime_error("No simulation has been loaded. Cannot construct simulation histogram.");
  initializeSimulationWeights(mainSimulation_,PionFluxWeighter_,KaonFluxWeighter_,PromptFluxWeighter_,osw_);
  simulation_initialized_=true;
}

bool Sterilizer::CheckSimulationInitialized() const {
  return(simulation_initialized_);
}

/*************************************************************************************************************
 * Functions to construct histograms
 * **********************************************************************************************************/

void Sterilizer::ConstructDataHistogram(){
  if(not data_loaded_))
    throw std::runtime_error("No data has been loaded. Cannot construct data histogram.");

  dataHist_ = HistType(LogarithmicAxis(0, 0.1), LinearAxis(0, 0.1), LinearAxis(2010, 1));

  dataHist_.getAxis(0)->setLowerLimit(minFitEnergy_);
  dataHist_.getAxis(0)->setUpperLimit(maxFitEnergy_);
  dataHist_.getAxis(1)->setLowerLimit(minCosth_);
  dataHist_.getAxis(1)->setUpperLimit(maxCosth_);

  // fill in the histogram with the data
  bin(sample_, dataHist_, binner);
  data_histogram_constructed_=true;
}

bool Sterilizer::CheckDataHistogramConstructed() const {
  return(data_histogram_constructed_);
}

void Sterilizer::ConstructSimulationHistogram(){
  if(not simulation_loaded_))
    throw std::runtime_error("No simulation has been loaded. Cannot construct simulation histogram.");
  if(not data_histogram_constructed_))
    throw std::runtime_error("Data histogram needs to be constructed before simulation histogram.");
  simHist_ = makeEmptyHistogramCopy(dataHist_);
  bin(mainSimulation_, simHist_, binner);
  simulation_histogram_constructed_=true;
}

bool Sterilizer::CheckSimulationHistogramConstructed() const {
  return(simulation_histogram_constructed_);
}

/*************************************************************************************************************
 * Functions to obtain distributions
 * **********************************************************************************************************/

marray<double,3> Sterilizer::GetDataDistribution() const {
    marray<double,3> array {static_cast<size_t>(dataHist_.getBinCount(2)),
                            static_cast<size_t>(dataHist_.getBinCount(1)),
                            static_cast<size_t>(dataHist_.getBinCount(0))};

    for(size_t iy=0; iy<data_hist.getBinCount(2); iy++){ // year
      for(size_t ic=0; ic<data_hist.getBinCount(1); ic++){ // zenith
        for(size_t ie=0; ie<data_hist.getBinCount(0); ie++){ // energy
          auto itc = static_cast<likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(dataHist_(ie,ic,iy));
          array[iy][ic][ie] = itc.size();
        }
      }
    }
    return array;
}


marray<double,3> Sterilizer::GetExpectation(SterileNeutrinoParameters snp, std::vector<double> nuisance) const {
    MakeSimulationHistogram(snp,nuisance);
    marray<double,3> array {static_cast<size_t>(simHist_.getBinCount(2)),
                            static_cast<size_t>(simHist_.getBinCount(1)),
                            static_cast<size_t>(simHist_.getBinCount(0))};

    auto weighter = DFWM(nuisance);
    for(size_t iy=0; iy<sim_hist.getBinCount(2); iy++){ // year
      for(size_t ic=0; ic<sim_hist.getBinCount(1); ic++){ // zenith
        for(size_t ie=0; ie<sim_hist.getBinCount(0); ie++){ // energy
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


marray<double,3> Sterilizer::GetRealization(SterileNeutrinoParameters snp, std::vector<double> nuisance) const{

}


// human interfaces
marray<double,3> Sterilizer::GetExpectation(SterileNeutrinoParameters snp, Nuisance nuisance)
{
  return GetExpectation(snp, ConvertNuisance(nuisance));
}

marray<double,3> Sterilizer::GetRealization(SterileNeutrinoParameters snp, Nuisance nuisance) const{
  return GetRealization(snp, ConvertNuisance(nuisance));
}

double Sterilizer::llhFull(SterileNeutrinoParameters snp, Nuisance nuisance) const {
  return llhFull(snp,ConvertNuisance(nuisance));
}




/*************************************************************************************************************
 * Functions to construct likelihood problem and evaluate it
 * **********************************************************************************************************/

void Sterilizer::ConstructLikelihoodProblem(Priors priors, Nuisance nuisanceSeed){
  if(not data_histogram_constructed_)
    throw std::runtime_error("Data histogram needs to be constructed before likelihood problem can be formulated.");
  if(not simulation_histogram_constructed_)
    throw std::runtime_error("Simulation histogram needs to be constructed before likelihood problem can be formulated.");
  
  llhpriors = ConvertPriors(priors);
  fitseed   = ConvertNuisance(nuisanceSeed);

  prob_ = likelihood::makeLikelihoodProblem<std::reference_wrapper<const Event>, 3, 6>(
      dataHist_, {simHist_}, llhpriors, {1.0}, likelihood::simpleDataWeighter(), DFWM,
      likelihood::poissonLikelihood(), fitseed );
  prob_.setEvaluationThreadCount(SteeringParams_.evalThreads);

  likelihood_problem_constructed_=true;
}

bool Sterilizer::CheckLikelihoodProblemConstruction() const {
return(likelihood_problem_constructed_);
}

double Sterilizer::EvalLLH(std::vector<double> nuisance) const {
  if(not likelihood_problem_constructed_)
    throw std::runtime_error("Likelihood problem has not been constructed..");
  return -prob_.evaluateLikelihood(nuisance);
}

<<<<<<< HEAD
fitResult Sterilizer::MinLLH(std::vector<std::pair<unsigned int,double>> fixedNuisanceParams) const {
=======
fitResult Sterilizer::MinLLH(Nuisance fixedParams) const {
>>>>>>> ed946e0c7263c123b5e64e7b572e6a46e02b1a00
  if(not likelihood_problem_constructed_)
    throw std::runtime_error("Likelihood problem has not been constructed..");
  
  std::vector<double> seed=prob.getSeed();
  std::vector<unsigned int> fixedIndices;
  
  std::vector<double> FixVec=ConvertNuisance(fixedParams);
  for(size_t i; i!=FixVec.size(); ++i)
      if(FixVec[i]>0.1) fixedIndices.push_back(i);
  

  return doFitLBFGSB(prob_, seed, fixedIndices);
}

/*************************************************************************************************************
 * Functions to change the sterile neutrino hypothesis
 * **********************************************************************************************************/

void Sterilizer::SetSterileNuParams(SterileNuParams snp){
  if(not simulation_loaded))
    throw std::runtime_error("No simulation has been loaded. Cannot weight to sterile hypothesis without simulation.");
  if(not mc_generation_weighter_constructed_)
    throw std::runtime_error("MonteCarlo generation weighter has to be constructed first.");
  if(not cross_section_weighter_constructed_)
    throw std::runtime_error("Cross section weighter has to be constructed first.");

  sterileNuParams_=snp;
  ConstructFluxWeighter();
  ConstructLeptonWeighter();
  WeightMC();
}

/*************************************************************************************************************
 * Functions to set options in the class
 * **********************************************************************************************************/
//Set the RNG seed
void Sterilizer::SetRandomNumberGeneratorSeed(unsigned int seed)
{
 steeringParams_.rngSeed=seed;
 rng_.seed(seed);
 if(!steeringParams_.quiet) std::cout<<"setting RNG seed to " << seed<<std::endl;
}


// Check that the directories where files are mean to be exist
bool Sterilizer::CheckDataPaths(DataPaths dp) const
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
CPrior Sterilizer::ConvertPriorSet(Priors pr)
{
  // construct continuous nuisance priors      
  UniformPrior  positivePrior(0.0,std::numeric_limits<double>::infinity());
  GaussianPrior normalizationPrior(pr.normCenter,pr.normWidth);
  GaussianPrior crSlopePrior(pr.crSlopePrior,pr.crSlopeWidth);
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

  if( delta_alpha.find(model_name) == delta_alpha.end() )
    throw std::runtime_error("Jordi delta key not found. Aborting.");
  double alpha = delta_alpha[delta_model_name];

  GaussianPrior ZCPrior(0.0,pr.zenithCorrectionMultiplier*alpha);

  // make and return priorset  
  return priors=makePriorSet(normalizationPrior,
                             positivePrior, 
                             positivePrior,
                             crSlopePrior,
                             simple_domEffPrior,
                             kaonPrior,
                             nanPrior,
                             ZCPrior);
}


// Given a human readable nuisance parameter set, make a nuisance vector
std::vector<double> ConvertNuisance(Nuisance ns)
{
  return std::vector<double> Nuisance{
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
