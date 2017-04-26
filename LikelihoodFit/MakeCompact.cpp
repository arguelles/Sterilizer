#import "SterileSearch.h"

int main(int argc, char* argv[]){
  using namespace SterileSearch;

  std::cout<<"====Making a compact data file ===="<<std::endl<<std::endl;

  // Initialize parameter sets at defaults
  DataPaths       dp;
  SteeringParams  sp;
  SterileNuParams snp;

  dp.squids_files_path="/data/user/bjones/Sterilizer/Sterilizer/conventional_fluxes/";
  dp.prompt_squids_files_path="/data/user/bjones/Sterilizer/Sterilizer/prompt_fluxes/";
  dp.compact_file_path="./";

  std::cout<<"==Making an object with non-compact data =="<<std::endl<<std::endl;
  sp.ReadCompact=false;
  std::shared_ptr<Sterilizer> ster = std::make_shared<Sterilizer>(dp, sp, snp);
  ster->ReportStatus();

  std::cout<<"Writing compact data"<<std::endl;
  ster->WriteCompact();

}

