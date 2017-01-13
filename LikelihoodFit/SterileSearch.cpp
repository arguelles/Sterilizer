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

void Sterilizer::LoadDOMEfficiencySplines(std::vector<int> years){
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

void Sterilizer::ConstructFluxWeighter(std::string conv_squids_files_path,std::string prompt_squids_files_path,std::string splines_path,SterileNeutrinoParameters snp){
  std::string sterile_neutrino_model_identifier = GetSterileNeutrinoModelIdentifier(snp);

  if(use_factorization_technique){
    std::string oscillation_model = "noint_atmospheric_"+sterile_neutrino_model_identifier+".hdf5";
    std::string atmospheric_model_pion = model_name + "_pion";
    std::string atmospheric_model_kaon = model_name + "_kaon";

    fluxPion_ = std::make_shared<LW::FactorizedSQUIDSFlux>(squids_files_path+oscillation_model,
                                                           flux_splines_path+atmospheric_model_pion+"_neutrino_spline.fits",
                                                           flux_splines_path+atmospheric_model_pion+"_antineutrino_spline.fits");
    fluxKaon_ = std::make_shared<LW::FactorizedSQUIDSFlux>(squids_files_path+oscillation_model,
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
  flux_prompt_ = std::make_shared<LW::SQUIDSFlux>(prompt_squids_files_path + "prompt_atmospheric_0.000000_0.000000.hdf5");
  flux_weighter_constructed_=true;
}

bool Sterilizer::CheckFluxWeighterConstructed() const {
  return(flux_weighter_constructed_);
}

void Sterilizer::ConstructMonteCarloGenerationWeighter(std::vector<std::string> simSetsToLoad){
  for( std::string sim_name : simSetsToLoad )
    mcw.addGenerationSpectrum(simInfo.find(sim_name)->second.details);
  mc_generation_weighter_constructed_=true;
}

bool Sterilizer::CheckFluxWeighterConstructed() const {
  return(flux_weighter_constructed_);
}

void Sterilizer::ConstructLeptonWeighter(){
  if(not mc_generation_weighter_constructed)
    throw std::runtime_error("MonteCarlo generation weighter has to be constructed first.");
  if(not flux_weighter_constructed)
    throw std::runtime_error("Flux weighter has to be constructed first.");
  if(not cross_section_weighter_constructed)
    throw std::runtime_error("Cross section weighter has to be constructed first.");

  PionFluxWeighter_ = LW::LeptonWeighter(flux_pion_,xsw_,mcw_);
  KaonFluxWeighter_ = LW::LeptonWeighter(flux_kaon_,xsw_,mcw_);
  PromptFluxWeighter_ = LW::LeptonWeighter(flux_promp_t,xsw_,mcw_);
  lepton_weighter_constructed_=true;
}

bool SterileSearch::CheckLeptonWeighterConstructed() const {
  return(lepton_weighter_constructed_);
}

/*************************************************************************************************************
 * Functions to initialize the MC weights
 * **********************************************************************************************************/

void Sterilizer::WeightMC(){
  if(not lepton_weighter_constructed)
    throw std::runtime_error("LeptonWeighter has to be constructed first.");
  if(not simulation_loaded))
    throw std::runtime_error("No simulation has been loaded. Cannot construct simulation histogram.");
  initializeSimulationWeights(mainSimulation_,PionFluxWeighter_,KaonFluxWeighter_,PromptFluxWeighter_,osw_dc_);
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
    marray<double,3> array {static_cast<size_t>(data_hist.getBinCount(2)),
                            static_cast<size_t>(data_hist.getBinCount(1)),
                            static_cast<size_t>(data_hist.getBinCount(0))};

    for(size_t iy=0; iy<data_hist.getBinCount(2); iy++){ // year
      for(size_t ic=0; ic<data_hist.getBinCount(1); ic++){ // zenith
        for(size_t ie=0; ie<data_hist.getBinCount(0); ie++){ // energy
          auto itc = static_cast<likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(data_hist(ie,ic,iy));
          array[iy][ic][ie] = itc.size();
        }
      }
    }
    return array;
}

marray<double,3> Sterilizer::GetExpectation(SterileNeutrinoParameters snp, std::vector<double> nuisance) const {
    MakeSimulationHistogram(snp,nuisance);
    marray<double,3> array {static_cast<size_t>(sim_hist.getBinCount(2)),
                            static_cast<size_t>(sim_hist.getBinCount(1)),
                            static_cast<size_t>(sim_hist.getBinCount(0))};

    auto weighter = DFWM(nuisance);
    for(size_t iy=0; iy<sim_hist.getBinCount(2); iy++){ // year
      for(size_t ic=0; ic<sim_hist.getBinCount(1); ic++){ // zenith
        for(size_t ie=0; ie<sim_hist.getBinCount(0); ie++){ // energy
          auto itc = static_cast<likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(sim_hist(ie,ic,iy));
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

/*************************************************************************************************************
 * Functions to construct histograms
 * **********************************************************************************************************/

double Sterilizer::llhFull(SterileNeutrinoParameters snp, std::vector<double> nuisance) const {

}

fitResult Sterilizer::llh(SterileNeutrinoParameters snp) const {

}




