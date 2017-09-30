#include "SterileSearch.h"

int main(int argc, char* argv[]){
  using namespace SterileSearch;

  // Initialize parameter sets at defaults
  DataPaths       dp;
  SteeringParams  sp;
  SterileNuParams snp;


  dp.squids_files_path="/data/user/bjones/Sterilizer/Sterilizer/test_data/";
  sp.readCompact=false;

  //Make the object
  std::shared_ptr<Sterilizer> ster = std::make_shared<Sterilizer>(dp, sp, snp);
 }
