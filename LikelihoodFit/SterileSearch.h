#ifndef _H_STERILE_SEARCH_
#define _H_STERILE_SEARCH_

struct fitResult {
  std::vector<double> params;
  double likelihood;
  unsigned int nEval, nGrad;
  bool succeeded;
};

#include <sys/types.h>
#include <sys/stat.h>



struct DataPaths {
  std::string compact_file_path =        "../compact_data/";
  std::string squids_files_path =        "/home/carguelles/work/TheSterileSearch/the_new_fluxes/";
  std::string prompt_squids_files_path = "/data/ana/NuFSGenMC/IC86_hypotheses/fluxes/prompt_0/";
  std::string xs_spline_path =           "/data/ana/NuFSGenMC/CrossSections/";
  std::string data_path =                "/data/ana/NuFSGenMC/Data/";
  std::string mc_path =                  "/data/ana/NuFSGenMC/Merged/";
  std::string oversize_function_main =   "/data/ana/NuFSGenMC/OversizeCorrections/NullCorrection.dat";
  std::string oversize_function_dc =     "/data/ana/NuFSGenMC/OversizeCorrections/NullCorrection.dat";
  std::string domeff_spline_path =       "/data/ana/NuFSGenMC/DomEffSplines/";
  std::string flux_splines_path =        "/home/carguelles/programs/SNOT/FluxSplines/propagated_splines/";
};

struct SteeringParams {
  double minFitEnergy=4e2;
  double maxFitEnergy=2.e4;
  double minCosth = -1.;
  double maxCosth = 0.2;
  uint32_t rngSeed=0;
  bool useFactorization=false;
  bool useBurnSample=false;
  std::string simToLoad="nufsgen_mie_0_99";
  std::string secondarySimToLoad="";
  bool quiet=false;
  std::string model_name = "HondaGaisser";
  std::string delta_model_name = "HondaGaisser";
  std::string oversizeFunction="NullCorrection";
  bool reuseBestFit=true;
  bool ReadCompact=true;
}

struct SterileNuParams {
  double th14 = 0;
  double th24 = 0;
  double th34 = 0;
  double del14 = 0;
  double del24 = 0;
  double dm41sq = 0;
  SterileNuParams(double th14,double th24, double th34, double del14, double del24, double dm41sq):
    th14(th14),th24(th24),th34(th34),del14(del14),del24(del24)
};

using namespace phys_tools::histograms;
using namespace likelihood;
using HistType = histogram<3,entryStoringBin<std::reference_wrapper<const Event>>>;

auto binner = [](HistType& h, const Event& e){
  h.add(e.energy,cos(e.zenith),e.year,amount(std::cref(e)));
};



class Sterilizer {
  private:
  // All the local configuration variables
    SteeringParams  steeringParams_;
    DataPaths       dataPaths_;

    //hypothesis point we fit to
    SterileNuParams sterileNuParams_;
    
    // To store best fit point, fit seed, and data challenge
    std::vector<double> existingBest_;
    std::vector<double> fitSeed_{1.02,0,0,0.05,.0985,1.1,1,0};
    std::vector<double> dataChallenge_nuisance_parameters_{1,0,0,0,.1,1,1,0};

    std::vector<double> livetime_;

    // to store events
    std::deque<Event> mainSimulation_;
    std::deque<Event> alternativeSimulation_;
    std::deque<Event> sample_;

    // random number generator
    std::mt19937 rng_;

    // histograms
    HistType dataHist_;
    HistType simHist_;

    // weighter object
    DiffuseFitWeighterMaker DFWM;
    LW::LeptonWeighter pionFluxWeighter_;
    LW::LeptonWeighter kaonFluxWeighter_;
    LW::LeptonWeighter promptFluxWeighter_;
  public:
    /// \brief Constructor
    Sterilizer(DataPaths dataPaths, SteeringParams steeringParams, SterileNeutrinoParameters snp);


    // Check that the directories where files are mean to be exist
    bool CheckDataPaths(DataPaths dp);

    // Check a directory exists and throw a relevant error otherwise.
    bool CheckDataPath(std::string p);
      
  protected:
    void LoadData(std::string filepath);
    void LoadCompactData(std::string filepath);
    void LoadMC(std::string filepath) {}
    void WeightMC(SterileNeutrinoParameters snp, std::vector<double> nuisance){}
    void LoadCompactMC(std::string filepath) {}
    void LoadFluxes(std::string filepath,SterileNeutrinoParameters snp) {}
    void MakeDataHistogram() {}
    void MakeSimulationHistogram(SterileNeutrinoParameters snp, std::vector<double> nuisance) {}

    bool CheckDataPaths(DataPaths dp);
    bool CheckDataPath(std::string p);

    void SetRandomNumberGeneratorSeed(unsigned int seed);




  public:
    marray<double,3> GetDataDistribution();
    marray<double,3> GetExpectation(SterileNeutrinoParameters snp, std::vector<double> nuisance);
    marray<double,3> GetRealization(SterileNeutrinoParameters snp, std::vector<double> nuisance);
    double llhFull(SterileNeutrinoParameters snp, std::vector<double> nuisance){}
    fitResult llh(SterileNeutrinoParameters snp) {}
    // set functions


    SteeringParams       GetSteeringParams()  { return steeringParams_; };
    DataPaths            GetDataPaths()       { return dataPaths_; };
    SterileNuParams      GetSterileNuParams() { return sterileNuParams_;};

    void       SetSteeringParams(SteeringParams p)   { steeringParams_=p; SetRandomNumberGeneratorSeed(p.rngSeed);};
    void       SetDataPaths(DataPaths p)             { dataPaths_=p; CheckDataPaths(p); };
    void       SetSterileNuParams(SterileNuParams p) { sterileNuParams_=p;};
};


Sterilizer::Sterilizer(DataPaths dataPaths, SteeringParams steeringParams, SterileNeutrinoParameters snp){
      // Set the parameter sets of this object
      steeringParams_=steeringParams;
      dataPaths_=dataPaths;
    
      // Set up the RNG
      SetRandomNumberGeneratorSeed(steeringParams.rngSeed);

      if(readCompact_){
        LoadCompactData(dataPaths_.compact_file_path);
        LoadCompactMC(dataPaths_.compact_file_path);
      } else {
        LoadData(dataPaths_.data_path);
        LoadMC(dataPaths_.mc_path);
      }
      LoadFluxes(dataPaths_.flux_splines_path,snp)
    }


void Sterilizer::SetRandomNumberGeneratorSeed(unsigned int seed)
{ 
  steeringParams_.rngSeed=seed;
  rng_.seed(seed);
  if(!steeringParams_.quiet) std::cout<<"setting RNG seed to " << seed<<std::endl;
}

// Check that the directories where files are mean to be exist
bool Sterilizer::CheckDataPaths(DataPaths dp)
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
bool Sterilizer::CheckDataPath(std::string p)
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

#endif
