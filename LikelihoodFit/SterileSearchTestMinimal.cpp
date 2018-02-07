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

  Nuisance n;
  auto expect = ster->GetExpectation(n);
  std::cout << expect.extent(0) << " " << expect.extent(1) << " " << expect.extent(2) << std::endl;

  for(size_t iy=0; iy<expect.extent(0); iy++){ // year
    for(size_t ic=0; ic<expect.extent(1); ic++){ // zenith
      for(size_t ie=0; ie<expect.extent(2); ie++){ // energy
        std::cout << iy << " "<< ic << " " << ie << " " << expect[iy][ic][ie] << std::endl;
      }
    }
  }

  auto result = ster->MinLLH();

  std::cout << result.likelihood << std::endl;

  std::cout << result.params.normalization << std::endl;
  /*
  auto expect_it = expect.begin();
  auto extents = expect.get_extents();
  for(auto it=extents.begin(); it<extents.end(); it++) {
      std::cout << *it << ' ';
  }
  std::cout << std::endl;
  */

 }
