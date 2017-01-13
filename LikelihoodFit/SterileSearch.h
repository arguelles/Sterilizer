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
#include "oversizeWeight.h"


struct DataPaths {
  std::string compact_file_path =        "../compact_data/";
  std::string squids_files_path =        "/home/carguelles/work/TheSterileSearch/the_new_fluxes/";
  std::string prompt_squids_files_path = "/data/ana/NuFSGenMC/IC86_hypotheses/fluxes/prompt_0/";
  std::string xs_spline_path =           "/data/ana/NuFSGenMC/CrossSections/";
  std::string data_path =                "/data/ana/NuFSGenMC/Data/";
  std::string mc_path =                  "/data/ana/NuFSGenMC/Merged/";
  std::string oversize_path =            "/data/ana/NuFSGenMC/OversizeCorrections/";
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
  std::string modelName = "HondaGaisser";
  std::string deltaModelName = "HondaGaisser";
  std::string oversizeFunction="NullCorrection";
  bool ReadCompact=true;
  std::string xs_model_name="";
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
    std::shared_ptr<LW::Flux> fluxKaon_,fluxPion_,fluxPrompt_;
    std::shared_ptr<LW::CrossSectionFromSpline> xsw_;
    LW::mcgenWeighter mcw_;
    OversizeWeighter osw_;

    bool cross_section_weighter_constructed_(false);
    bool flux_section_weighter_constructed_(false);
    bool lepton_weighter_constructed_(false);
    bool oversize_weighter_constructed_(false);

    // DOM efficiency splines
    std::vector<std::unique_ptr<Splinetable>> domEffConv_;
    std::vector<std::unique_ptr<Splinetable>> domEffPrompt_;

  public:
    // Constructor
    Sterilizer(DataPaths dataPaths, SteeringParams steeringParams, SterileNeutrinoParameters snp);

    // Check that the directories where files are mean to be exist
    bool CheckDataPaths(DataPaths dp);

    // Check a directory exists and throw a relevant error otherwise.
    bool CheckDataPath(std::string p);
      
  protected:
    // Functions to load and unload data
    void LoadData(std::string filepath);
    void LoadMC(std::string filepath) {}
    void LoadCompact(std::string filepath);
    void WriteCompact(std::string filepath);
    // Functions to load and unload data
    void WeightMC(SterileNeutrinoParameters snp, std::vector<double> nuisance){}
    void LoadFluxes(std::string filepath,SterileNeutrinoParameters snp) {}
    void MakeDataHistogram() {}
    void MakeSimulationHistogram(SterileNeutrinoParameters snp, std::vector<double> nuisance) {}

    bool CheckDataPaths(DataPaths dp);
    bool CheckDataPath(std::string p);

    void SetRandomNumberGeneratorSeed(unsigned int seed);

 private:
    // Test if two particles are in the same generation (e.g., mu / numu)
    bool SameGeneration(particleType p1, particleType p2);

    // Make the fit
    template<typename LikelihoodType> fitResult DoFitLBFGSB(LikelihoodType& likelihood, const std::vector<double>& seed,std::vector<unsigned int> indicesToFix);

    void ConstructFluxWeighter(std::string squids_files_path,std::string splines_path,SterileNeutrinoParameters snp);

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





#endif
