#import "SterileSearch.h"

int main(int argc, char* argv[]){
  using namespace SterileSearch;

  std::cout<<"==== Running the sterilizer test module ===="<<std::endl<<std::endl;

  // Initialize parameter sets at defaults
  DataPaths       dp;
  SteeringParams  sp;
  SterileNuParams snp;

  //dp.mc_path="/data/user/bjones/Sterilizer/Sterilizer/test_data/";
  //dp.squids_files_path="/data/user/bjones/Sterilizer/Sterilizer/test_data/";

  dp.mc_path="../mc/";
  dp.squids_files_path="../conventional_fluxes/";
  dp.prompt_squids_files_path="../prompt_fluxes/";
  dp.domeff_spline_path="../DOMEfficiencyFits/";
  dp.xs_spline_path="../cross_sections/";
  dp.data_path="../data/";
  dp.oversize_path="../oversize_corrections/";
  dp.flux_splines_path="../propagated_splines/";

  std::cout<<"==Making an object with non-compact data =="<<std::endl<<std::endl;
  sp.ReadCompact=false;

  //Make the object
  std::shared_ptr<Sterilizer> ster = std::make_shared<Sterilizer>(dp, sp, snp);
  ster->ReportStatus();

  std::cout<< "Getting the data distribution"<<std::endl;
  marray<double,3> DataDist=ster->GetDataDistribution();
  std::cout<<"   Bin 1,1 " << DataDist[0][1][1] << std::endl;
  std::cout<<"   Bin 5,5 " << DataDist[0][5][5] << std::endl;

  std::cout<< "Getting expectation for nominal nuisance"<<std::endl;
  Nuisance ns;
  std::cout<< " Normalization : " <<ns.normalization<<std::endl;

  marray<double,3> NominalExpectation=ster->GetExpectation(ns);
  std::cout<<"   Bin 1,1 " << NominalExpectation[0][1][1] << std::endl;
  std::cout<<"   Bin 5,5 " << NominalExpectation[0][5][5] << std::endl;

  std::cout<< "Getting expectation for adjusted nuisance"<<std::endl;
  ns.normalization=1.2;
  std::cout<< " Normalization : " <<ns.normalization<<std::endl;
  marray<double,3> AdjustedExpectation=ster->GetExpectation(ns);
  std::cout<<"   Bin 1,1 " << AdjustedExpectation[0][1][1] << std::endl;
  std::cout<<"   Bin 5,5 " << AdjustedExpectation[0][5][5] << std::endl;

  std::cout<<"Getting a realization for adjusted nuisance"<<std::endl;
  int seed=100;
  marray<double,3> Realization=ster->GetRealization(ns,seed);
  /* // are we actiually filling it in.
  std::cout << Realization.extent(0) << std::endl;
  std::cout << Realization.extent(1) << std::endl;
  std::cout << Realization.extent(2) << std::endl;
  */
  std::cout<<"   Bin 1,1 " << Realization[0][1][1] << std::endl;
  std::cout<<"   Bin 5,5 " << Realization[0][5][5] << std::endl;

  ns.normalization=1.0;
  std::cout<< " Normalization : " <<ns.normalization<<std::endl;

  std::cout<<"Writing compact data"<<std::endl;
  ster->WriteCompact();

  std::cout<<std::endl<<"==Making an object with compact data=="<<std::endl<<std::endl;
  sp.ReadCompact=true;

  // share pointer takes care of deleting the previous one. Just go.
  ster = std::make_shared<Sterilizer>(dp, sp, snp);
  ster->ReportStatus();

  std::cout<< "Getting the data distribution"<<std::endl;
  DataDist=ster->GetDataDistribution();
  std::cout<<"   Bin 1,1 " << DataDist[0][1][1] << std::endl;
  std::cout<<"   Bin 5,5 " << DataDist[0][5][5] << std::endl;

  std::cout<<"Adjusting the sterile params to (1.0, 0.2)"<<std::endl;
  snp.th24=0.2;
  snp.dm41sq=1.0;
  ster->SetSterileNuParams(snp);

  std::cout<< "Getting the data distribution"<<std::endl;
  DataDist=ster->GetDataDistribution();
  std::cout<<"   Bin 1,1 " << DataDist[0][1][1] << std::endl;
  std::cout<<"   Bin 5,5 " << DataDist[0][5][5] << std::endl;

  std::cout<< "Getting expectation for nominal nuisance"<<std::endl;
  NominalExpectation=ster->GetExpectation(ns);
  std::cout<<"   Bin 1,1 " << NominalExpectation[0][1][1] << std::endl;
  std::cout<<"   Bin 5,5 " << NominalExpectation[0][5][5] << std::endl;

  std::cout<<std::endl<<"==Making an object with non-null sterile parameters from scratch=="<<std::endl<<std::endl;
  Sterilizer * stersig = new Sterilizer(dp, sp, snp);
  stersig->ReportStatus();

  std::cout<< "Getting the data distribution"<<std::endl;
  DataDist=stersig->GetDataDistribution();
  std::cout<<"   Bin 1,1 " << DataDist[0][1][1] << std::endl;
  std::cout<<"   Bin 5,5 " << DataDist[0][5][5] << std::endl;

  std::cout<< "Getting expectation for nominal nuisance"<<std::endl;
  NominalExpectation=ster->GetExpectation(ns);
  std::cout<<"   Bin 1,1 " << NominalExpectation[0][1][1] << std::endl;
  std::cout<<"   Bin 5,5 " << NominalExpectation[0][5][5] << std::endl;

  std::cout<<"==== Sterilizer test module complete ===="<<std::endl<<std::endl;
}

