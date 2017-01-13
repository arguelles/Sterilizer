#ifndef _H_STERILE_SEARCH_
#define _H_STERILE_SEARCH_

#include <sys/types.h>
#include <sys/stat.h>
#include <iterator>
#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <chrono>
#include <queue>
#include <vector>

#include "oversizeWeight.h"
#include "likelihood.h"
#include "Event.h"
#include "analysisWeighting.h"
#include "dataIO.h"
#include "runspec.h"
#include "oversizeWeight.h"


struct Nuisance {
  float normaliztion;
  float astroFlux;
  float promptFlux;
  float crSlope;
  float domEfficiency;
  float piKRatio;
  float nuNubarRatio;
  float zenithCorrection;
};

struct NuisanceFlag {
  bool normaliztion=0;
  bool astroFlux=1;
  bool promptFlux=1;
  bool crSlope=0;
  bool domEfficiency=0;
  bool piKRatio=0;
  bool nuNubarRatio=0;
  bool zenithCorrection=0;
};

struct Priors {
  float normCenter=1.;
  float normWidth=0.4;
  float crSlopeCenter=0.0;
  float crSlopeWidth=0.05;
  float domEffCenter=-0.1;
  float domEffWidth=0.3;
  float piKRatioCenter=1.0;
  float piKRatioWidth=0.1;
  float nuNubarRatioCenter=1.0;
  float nuNubarRatioWidth=0.1;
  float zenithCorrectionMultiplier=0.038;
};

struct FitResult {
  Nuisance params;
  double likelihood;
  unsigned int nEval, nGrad;
  bool succeeded;
};


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
  bool quiet=false;
  size_t evalThreads=1;
  std::string modelName = "HondaGaisser";
  std::string oversizeFunction="NullCorrection";
  bool ReadCompact=true;
  std::string xs_model_name="";
};

struct SterileNuParams {
  unsigned int modelId =0;
  double th14 = 0;
  double th24 = 0;
  double th34 = 0;
  double del14 = 0;
  double del24 = 0;
  double dm41sq = 0;
  SterileNuParams(unsigned int modelId, double th14,double th24, double th34, double del14, double del24, double dm41sq):
    th14(th14),th24(th24),th34(th34),del14(del14),del24(del24){}
};

using namespace nusquids;
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

    // Status flags
    bool cross_section_weighter_constructed_ = (false);
    bool flux_section_weighter_constructed_ = (false);
    bool lepton_weighter_constructed_ = (false);
    bool oversize_weighter_constructed_ = (false);
    bool dom_efficiency_splines_loaded_ = (false);
    bool data_histogram_constructed_ = (false);
    bool simulation_loaded_ = (false);
    bool mc_generation_weighter_constructed_ = (false);
    bool data_loaded_ = (false);

    // DOM efficiency splines
    std::vector<std::unique_ptr<Splinetable>> domEffConv_;
    std::vector<std::unique_ptr<Splinetable>> domEffPrompt_;

  public:
    // Constructor
    Sterilizer(DataPaths dataPaths, SteeringParams steeringParams, SterileNuParams snp);

    // Check that the directories where files are mean to be exist
    bool CheckDataPaths(DataPaths dp) const;

    // Check a directory exists and throw a relevant error otherwise.
    bool CheckDataPath(std::string p) const;

  protected:
    // Functions to load and unload data
    void LoadData();
    void LoadMC();
    void LoadCompact();
    void WriteCompact() const;
    void ClearData();
    void ClearSimulation();
    // loading DOM efficiency splines
    void LoadDOMEfficiencySplines(std::vector<unsigned int> years);
    // functions to construct the weighters
    void ConstructCrossSectionWeighter();
    void ConstructFluxWeighter();
    void ConstructMonteCarloGenerationWeighter();
    void ConstructLeptonWeighter();
    void ConstructOversizeWeighter();
    // Function to initialize the MC weights
    void WeightMC();
    // functions to construct the histograms of data and simulation
    void ConstructDataHistogram();
    void ConstructSimulationHistogram();
    // functions to construct the likelihood problem
    void ConstructLikelihoodProblem();

    // Converters between human and vector forms
    CPrior ConvertPriorSet(Priors pr);
    std::vector<double> ConvertNuisance(Nuisance ns);
    std::vector<bool> ConvertNuisanceFlag(NuisanceFlag ns);
    std::vector<double> ConvertVecToNuisance(std::vector<double> vecns);


  public:
    // functions to check the status of the object
    bool CheckDataLoaded() const;
    bool CheckSimulationLoaded() const;
    bool CheckDOMEfficiencySplinesLoaded() const;
    bool CheckCrossSectionWeighterConstructed() const;
    bool CheckFluxWeighterConstructed() const;
    bool CheckLeptonWeighterConstructed() const;
    bool CheckSimulationInitialized() const;
    bool CheckDataHistogramConstructed() const;
    bool CheckSimulationHistogramConstructed() const;
    bool CheckLikelihoodProblemConstruction() const;
    // functions to obtain distributions
    marray<double,3> GetDataDistribution() const;
    marray<double,3> GetExpectation(std::vector<double> nuisance) const;
    marray<double,3> GetExpectation(Nuisance nuisance) const;
    marray<double,3> GetRealization(std::vector<double> nuisance) const;
    marray<double,3> GetRealization(Nuisance nuisance) const;
    // functions to evaluate the likelihood
    double EvalLLH(std::vector<double> nuisance) const;
    double EvalLL(Nuisance nuisance) const;
    FitResult MinLLH(NuisanceFlag fixedParams) const;
    void SetSterileNuParams(SterileNuParams snp);
  private:
    // Make the fit
    template<typename LikelihoodType> FitResult DoFitLBFGSB(LikelihoodType& likelihood, const std::vector<double>& seed,std::vector<unsigned int> indicesToFix);
  public:
    // set functions

    SteeringParams       GetSteeringParams()  { return steeringParams_; };
    DataPaths            GetDataPaths()       { return dataPaths_; };
    SterileNuParams      GetSterileNuParams() { return sterileNuParams_;};

    void       SetSteeringParams(SteeringParams p)   { steeringParams_=p; SetRandomNumberGeneratorSeed(p.rngSeed);};
    void       SetDataPaths(DataPaths p)             { dataPaths_=p; CheckDataPaths(p); };
    void       SetRandomNumberGeneratorSeed(unsigned int seed);
};

#endif
