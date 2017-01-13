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

void SterileSearch::LoadData(std::string dataPath){
		try{
			sample=loadExperimentalData(dataPath,UseBurnsample);
		} catch(std::exception& ex){
			std::cerr << "Problem loading experimental data: " << ex.what() << std::endl;
			return(1);
		}
		if(!quiet)
			std::cout << "Loaded " << sample.size() << " experimental events" << std::endl;
}

void SterileSearch::LoadMC(std::string simulationPath,std::vector<std::string> simSetsToLoad){
    bool loadTargeted=true;
		try{
      loadSimulatedData(mainSimulation,simulationPath,livetime,simInfo,simSetsToLoad,loadTargeted);
		} catch(std::exception& ex){
			std::cerr << "Problem loading simulated data: " << ex.what() << std::endl;
			return(1);
		}
		if(!quiet)
			std::cout << "Loaded " << mainSimulation.size() << " events in main simulation set" << std::endl;
}

void SterileSearch::LoadCompact(std::string compact_data_path, std::string simulation_to_load) {
		try{
			unsplatData(compact_data_path+"/"+simulation_to_load+"_compact_data.dat",getFileChecksum(argv[0]),sample,mainSimulation);
			if(!quiet){
				std::cout << "Loaded " << sample.size() << " experimental events." << std::endl;
				std::cout << "Loaded " << mainSimulation.size() << " events in main simulation set." << std::endl;
			}
		} catch(std::runtime_error& re){
			std::cerr << re.what() << std::endl;
			std::cerr << "Failed to load compact data" << std::endl;
			return(1);
		}
}

void SterileSearch::WriteCompact(std::string compact_data_path, std::string simulation_to_load) const {
		try{
			splatData(compact_data_path+"/"+simulation_to_load+"_compact_data.dat",
                getFileChecksum(argv[0]),sample,mainSimulation);
		} catch(std::runtime_error& re){
			std::cerr << re.what() << std::endl;
			std::cerr << "Failed to save compact data" << std::endl;
			return(1);
		}
}

/*************************************************************************************************************
 * Functions to load nusquids fluxes
 * **********************************************************************************************************/

void SterileSearch::LoadFluxes(std::string filepath,SterileNeutrinoParameters snp) {
}

/*************************************************************************************************************
 * Functions to load to load DOM efficiency splines
 * **********************************************************************************************************/

void SterileSearch::LoadDOMEfficiencySplines(){
  for(size_t year_index=0; year_index<years.size(); year_index++){
    domEffConv[year_index] = std::unique_ptr<Splinetable>(new Splinetable(domeff_spline_path+"/conv_IC"+std::to_string(years[year_index]+".fits"));
    domEffPrompt[year_index] = std::unique_ptr<Splinetable>(new Splinetable(domeff_spline_path+"/prompt_IC"+std::to_string(years[year_index])+".fits"));
  }
}

/*************************************************************************************************************
 * Functions to construct weighters
 * **********************************************************************************************************/

void SterileSearch::ConstructCrossSectionWeighter(std::string xs_spline_path, std::string xs_model_name){
  xsw = std::make_shared<LW::CrossSectionFromSpline>(static_cast<std::string>(xs_spline_path),xs_model_name);
  cross_section_weighter_constructed=true;
}

void SterileSearch::ConstructFluxWeighter(std::string conv_squids_files_path,std::string prompt_squids_files_path,std::string splines_path,SterileNeutrinoParameters snp){
  std::string sterile_neutrino_model_identifier = GetSterileNeutrinoModelIdentifier(snp);

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
  flux_prompt = std::make_shared<LW::SQUIDSFlux>(prompt_squids_files_path + "prompt_atmospheric_0.000000_0.000000.hdf5");
  flux_weighter_constructed=true;
}

void SterileSearch::ConstructMonteCarloGenerationWeighter(std::vector<std::string> simSetsToLoad){
  for( std::string sim_name : simSetsToLoad )
    mcw.addGenerationSpectrum(simInfo.find(sim_name)->second.details);
  mc_generation_weighter_constructed=true;
}

void SterileSearch::ConstructLeptonWeighter(){
  if(not mc_generation_weighter_constructed)
    throw std::runtime_error("MonteCarlo generation weighter has to be constructed first.");
  if(not flux_weighter_constructed)
    throw std::runtime_error("Flux weighter has to be constructed first.");
  if(not cross_section_weighter_constructed)
    throw std::runtime_error("Cross section weighter has to be constructed first.");

  PionFluxWeighter = LW::LeptonWeighter(flux_pion,xsw,mcw);
  KaonFluxWeighter = LW::LeptonWeighter(flux_kaon,xsw,mcw);
  PromptFluxWeighter = LW::LeptonWeighter(flux_prompt,xsw,mcw);
  mc_generation_weighter_constructed=true;
}

/*************************************************************************************************************
 * Functions to obtain distributions
 * **********************************************************************************************************/

marray<double,3> SterileSearch::GetDataDistribution() const {
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

marray<double,3> SterileSearch::GetExpectation(SterileNeutrinoParameters snp, std::vector<double> nuisance) const {
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

marray<double,3> SterileSearch::GetRealization(SterileNeutrinoParameters snp, std::vector<double> nuisance) const{

}
