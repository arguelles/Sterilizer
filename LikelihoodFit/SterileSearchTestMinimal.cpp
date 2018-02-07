#include "SterileSearch.h"

int main(int argc, char* argv[]){
  using namespace SterileSearch;

  // Initialize parameter sets at defaults
  DataPaths       dp;
  SteeringParams  sp;
  SterileNuParams snp;

  dp.squids_files_path="/data/user/bjones/Sterilizer/Sterilizer/test_data/";
  dp.data_path="/data/ana/NuFSGenMC/Data/";
  dp.data_path="/data/ana/NuFSGenMC/Pass2_IC86.2011/Data/";
  dp.mc_path="/data/ana/NuFSGenMC/Merged/";
  dp.prompt_squids_files_path="/home/carguelles/work/Sterilizer/prompt_fluxes/";

  sp.readCompact=false;
  sp.calculate_nusquids_on_the_fly=true;
  sp.quiet=false;
  sp.modelName="PolyGonato_QGSJET-II-04";
  sp.use_simplified_simulation = true;

  //Make the object
  std::shared_ptr<Sterilizer> ster = std::make_shared<Sterilizer>(dp, sp, snp);
 }
