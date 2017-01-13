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
    sample=loadExperimentalData(dataPaths.data_path,steeringParams_.useBurnsample);
  } catch(std::exception& ex){
    std::cerr << "Problem loading experimental data: " << ex.what() << std::endl;
    return(1);
  }
  if(!quiet)
    std::cout << "Loaded " << sample.size() << " experimental events" << std::endl;
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
}

void Sterilizer::LoadCompact() {
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

/*************************************************************************************************************
 * Functions to load nusquids fluxes
 * **********************************************************************************************************/

void Sterilizer::LoadFluxes(std::string filepath,SterileNeutrinoParameters snp) {
}

/*************************************************************************************************************
 * Functions to load to load DOM efficiency splines
 * **********************************************************************************************************/

void Sterilizer::LoadDOMEfficiencySplines(std::vector<int> years){
  for(size_t year_index=0; year_index<years.size(); year_index++){
    domEffConv_[year_index] = std::unique_ptr<Splinetable>(new Splinetable(dataPaths_.domeff_spline_path+"/conv_IC"+std::to_string(years[year_index]+".fits"));
    domEffPrompt_[year_index] = std::unique_ptr<Splinetable>(new Splinetable(dataPaths_.domeff_spline_path+"/prompt_IC"+std::to_string(years[year_index])+".fits"));
  }
}

/*************************************************************************************************************
 * Functions to construct weighters
 * **********************************************************************************************************/

void Sterilizer::ConstructCrossSectionWeighter(){
  xsw_ = std::make_shared<LW::CrossSectionFromSpline>(static_cast<std::string>(dataPaths.xs_spline_path),steeringParams_.xs_model_name);
  cross_section_weighter_constructed_=true;
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
      fluxKaon = std::make_shared<LW::SQUIDSFlux>(dataPaths_.squids_files_path + flux_kaon_filename + ".hdf5");
      fluxPion = std::make_shared<LW::SQUIDSFlux>(dataPaths_.squids_files_path + flux_pion_filename + ".hdf5");
  }
  fluxPrompt = std::make_shared<LW::SQUIDSFlux>(dataPaths_.prompt_squids_files_path + "prompt_atmospheric_0.000000_0.000000.hdf5");
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
  if(not cross_section_weighter_constructed_)
    throw std::runtime_error("Cross section weighter has to be constructed first.");

  PionFluxWeighter_ = LW::LeptonWeighter(fluxPion_,xsw_,mcw_);
  KaonFluxWeighter_ = LW::LeptonWeighter(fluxKaon_,xsw_,mcw_);
  PromptFluxWeighter_ = LW::LeptonWeighter(fluxPrompt_,xsw_,mcw_);
  lepton_weighter_constructed_=true;
}

/*************************************************************************************************************
 * Functions to obtain distributions
 * **********************************************************************************************************/

void Sterilizer::WeightMC(){
  if(not lepton_weighter_constructed_)
    throw std::runtime_error("LeptonWeighter has to be constructed first.");

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
