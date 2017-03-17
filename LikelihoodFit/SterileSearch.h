#ifndef _H_STERILE_SEARCH_
#define _H_STERILE_SEARCH_
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <iterator>
#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <chrono>
#include <queue>
#include <vector>
#include <memory>

#include "likelihood.h"
#include "Event.h"
#include "analysisWeighting.h"
#include "compactIO.h"
#include "runspec.h"
#include "oversizeWeight.h"

namespace SterileSearch {

struct Nuisance {
  float normalization=1.0;
  float astroFlux=0.0;
  float promptFlux=0.0;
  float crSlope=0.0;
  float domEfficiency=0.1;
  float piKRatio=1.0;
  float nuNubarRatio=1.0;
  float zenithCorrection=0.0;
  Nuisance(){};
};

struct NuisanceFlag {
  bool normalization=false;
  bool astroFlux=true;
  bool promptFlux=true;
  bool crSlope=false;
  bool domEfficiency=false;
  bool piKRatio=false;
  bool nuNubarRatio=false;
  bool zenithCorrection=false;
  NuisanceFlag(){};
};

struct Priors {
  float normCenter=1.;
  float normWidth=0.4;
  float crSlopeCenter=0.0;
  float crSlopeWidth=0.05;
  float domEffCenter=0.1;
  float domEffWidth=0.1;
  float piKRatioCenter=1.0;
  float piKRatioWidth=0.1;
  float nuNubarRatioCenter=1.0;
  float nuNubarRatioWidth=0.1;
  float zenithCorrectionMultiplier=0.038;
  Priors(){};
};

struct FitResult {
  Nuisance params;
  double likelihood;
  unsigned int nEval, nGrad;
  bool succeeded;
  FitResult(){};
};

struct DataPaths {
  std::string compact_file_path =        "../compact_data/";
  std::string squids_files_path =        "/data/user/twatson/fluxes/flux_0/";
  std::string prompt_squids_files_path = "/data/ana/NuFSGenMC/IC86_hypotheses/fluxes/prompt_0/";
  std::string xs_spline_path =           "/data/ana/NuFSGenMC/CrossSections/";
  std::string data_path =                "/data/ana/NuFSGenMC/Data/";
  std::string mc_path =                  "/data/ana/NuFSGenMC/Merged/";
  std::string oversize_path =            "/data/ana/NuFSGenMC/OversizeCorrections/";
  std::string domeff_spline_path =       "/data/ana/NuFSGenMC/DomEffSplines/";
  std::string flux_splines_path =        "/home/carguelles/programs/SNOT/FluxSplines/propagated_splines/";
  DataPaths(){};
};


struct SteeringParams {
  float minFitEnergy=4e2;
  float maxFitEnergy=2.e4;
  float minCosth = -1.;
  float maxCosth = 0.2;
  float logEbinEdge = 2.6;
  float logEbinWidth = 0.169;
  float cosThbinEdge = 0.0;
  float cosThbinWidth = 0.06;
  bool useFactorization=false;
  bool useBurnSample=false;
  std::string simToLoad="nufsgen_mie_0_99";
  bool quiet=false;
  size_t evalThreads=1;
  std::string modelName = "HondaGaisser";
  std::string oversizeFunction="NullCorrection";
  bool ReadCompact=true;
  std::string xs_model_name="";
  std::vector<unsigned int> years={2011};
  std::map<unsigned int, double> burnSampleLivetime = std::map<unsigned int,double>{{2011,758.59*60*60}};
  std::map<unsigned int, double> fullLivetime= std::map<unsigned int,double>{{2011,8249.6*3600}};

  SteeringParams(){};
};

struct SterileNuParams {
  unsigned int modelId = 0;
  double th14 = 0;
  double th24 = 0;
  double th34 = 0;
  double del14 = 0;
  double del24 = 0;
  double dm41sq = 0;
  SterileNuParams(unsigned int modelId, double th14,double th24, double th34, double del14, double del24, double dm41sq):
    th14(th14),th24(th24),th34(th34),del14(del14),del24(del24){}
  SterileNuParams(){};
};

using namespace nusquids;
using namespace phys_tools::histograms;
using namespace likelihood;
using HistType = histogram<3,entryStoringBin<std::reference_wrapper<const Event>>>;
using CPrior=FixedSizePriorSet<GaussianPrior,UniformPrior,UniformPrior,GaussianPrior,LimitedGaussianPrior,GaussianPrior,GaussianPrior,GaussianPrior>;
using LType=LikelihoodProblem<std::reference_wrapper<const Event>,simpleDataWeighter,DiffuseFitWeighterMaker,CPrior,poissonLikelihood,3,8>;

template<typename ContainerType, typename HistType, typename BinnerType>
  void bin(const ContainerType& data, HistType& hist, const BinnerType& binner){
  for(const Event& event : data)
    binner(hist,event);
}

class Sterilizer {
  private:
  // All the local configuration variables
    SteeringParams  steeringParams_;
    DataPaths       dataPaths_;

    //hypothesis point we fit to
    SterileNuParams sterileNuParams_;

    // to store events
    std::deque<Event> mainSimulation_;
    std::deque<Event> sample_;

    // histograms
    HistType dataHist_;
    HistType simHist_;

    // minimizing objects
    Nuisance     fitSeed_;
    NuisanceFlag fixedParams_;

    // weighter object
    DiffuseFitWeighterMaker DFWM;
    std::shared_ptr<LW::LeptonWeighter> pionFluxWeighter_;
    std::shared_ptr<LW::LeptonWeighter> kaonFluxWeighter_;
    std::shared_ptr<LW::LeptonWeighter> promptFluxWeighter_;
    std::shared_ptr<LW::Flux> fluxKaon_,fluxPion_,fluxPrompt_;
    std::shared_ptr<LW::CrossSectionFromSpline> xsw_;
    LW::mcgenWeighter mcw_;
    OversizeWeighter osw_;

    // Status flags
    bool xs_weighter_constructed_ = (false);
    bool flux_weighter_constructed_ = (false);
    bool lepton_weighter_constructed_ = (false);
    bool oversize_weighter_constructed_ = (false);
    bool dom_efficiency_splines_constructed_ = (false);
    bool data_histogram_constructed_ = (false);
    bool simulation_histogram_constructed_ = (false);
    bool simulation_loaded_ = (false);
    bool mc_generation_weighter_constructed_ = (false);
    bool data_loaded_ = (false);
    bool likelihood_problem_constructed_ = (false);
    bool simulation_initialized_ = (false);

    // DOM efficiency splines
    std::map<unsigned int,std::shared_ptr<Splinetable>> domEffConv_;
    std::map<unsigned int,std::shared_ptr<Splinetable>> domEffPrompt_;

    // likehood problem object
    std::shared_ptr<LType> prob_;
  public:
    // Constructor
    Sterilizer(DataPaths dataPaths, SteeringParams steeringParams, SterileNuParams snp);

    // Check that the directories where files are mean to be exist
    bool CheckDataPaths(DataPaths dp) const;

    // Check a directory exists and throw a relevant error otherwise.
    bool CheckDataPath(std::string p) const;

    // Error trap bad file paths
    std::string CheckedFilePath(std::string) const;

    void WriteCompact() const;

  protected:
    // Functions to load and unload data
    void LoadData();
    void LoadMC();
    void LoadCompact();
    void ClearData();
    void ClearSimulation();
    // loading DOM efficiency splines
    void LoadDOMEfficiencySplines();
    // functions to construct the weighters
    void ConstructCrossSectionWeighter();
    void ConstructFluxWeighter();
    void ConstructMonteCarloGenerationWeighter();
    void ConstructLeptonWeighter();
    void ConstructOversizeWeighter();
    // Function to initialize the MC weights
    void WeightMC();
    void InitializeSimulationWeights();
    // functions to construct the histograms of data and simulation
    void ConstructDataHistogram();
    void ConstructSimulationHistogram();
    // functions to construct the likelihood problem
    void ConstructLikelihoodProblem(Priors priors, Nuisance nuisanceSeed, NuisanceFlag fixedParams);

    double GetZenithCorrectionScale() const;
    // Converters between human and vector forms
    std::vector<double> ConvertNuisance(Nuisance ns) const;
    std::vector<bool> ConvertNuisanceFlag(NuisanceFlag ns) const;
    Nuisance ConvertVecToNuisance(std::vector<double> vecns) const;
    marray<double,3> GetRealization(std::vector<double> nuisance, int seed) const;
    marray<double,3> GetExpectation(std::vector<double> nuisance) const;
    marray<double,2> SpitRealization(std::vector<double> nuisance, int seed) const;
    marray<double,2> SpitExpectation( std::vector<double> nuisance) const;
    std::vector<double> PullBinEdges(int dim, const  HistType& h) const;
    void SetupAsimov(std::vector<double> Nuisance);

  public:

    // Methods to spit out and swallow event samples
    marray<double,2> SpitData() const;
    marray<double,2> SpitRealization(Nuisance nuisance, int seed) const;
    marray<double,2> SpitExpectation(Nuisance nuisance) const;
    double Swallow(marray<double,2> Data);
    bool SetupAsimov(Nuisance nuisance);

    // Methods to get histogram binning
    std::vector<double> GetEnergyBinsData() const;
    std::vector<double> GetZenithBinsData() const;
    std::vector<double> GetEnergyBinsMC() const;
    std::vector<double> GetZenithBinsMC() const;

    // functions to check the status of the object
    bool CheckDataLoaded() const                       {return data_loaded_;};
    bool CheckSimulationLoaded() const                 {return simulation_loaded_;};
    bool CheckDOMEfficiencySplinesConstructed() const  {return dom_efficiency_splines_constructed_;};
    bool CheckCrossSectionWeighterConstructed() const  {return xs_weighter_constructed_;};
    bool CheckFluxWeighterConstructed() const          {return flux_weighter_constructed_;};
    bool CheckOversizeWeighterConstructed() const      {return oversize_weighter_constructed_;};
    bool CheckLeptonWeighterConstructed() const        {return lepton_weighter_constructed_;};
    bool CheckDataHistogramConstructed() const         {return data_histogram_constructed_;};
    bool CheckSimulationHistogramConstructed() const   {return simulation_histogram_constructed_;};
    bool CheckLikelihoodProblemConstruction() const    {return likelihood_problem_constructed_;};
    void ReportStatus() const;

    // functions to obtain distributions
    marray<double,3> GetDataDistribution() const;
    marray<double,3> GetExpectation(Nuisance nuisance) const;
    marray<double,3> GetRealization(Nuisance nuisance, int seed) const;
    // functions to evaluate the likelihood
    double EvalLLH(std::vector<double> nuisance) const;
    double EvalLLH(Nuisance nuisance) const;
    FitResult MinLLH() const;
    void SetSterileNuParams(SterileNuParams snp);

  private:

    // Nasty template part of fit function
    template<typename LikelihoodType>
      bool DoFitLBFGSB(LikelihoodType& likelihood, LBFGSB_Driver& minimizer) const{
      using namespace likelihood;
      bool succeeded=
	minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));
      return succeeded;
    };

  public:
    // set functions

    SteeringParams       GetSteeringParams()  { return steeringParams_; };
    DataPaths            GetDataPaths()       { return dataPaths_; };
    SterileNuParams      GetSterileNuParams() { return sterileNuParams_;};

    void       SetSteeringParams(SteeringParams p)   { steeringParams_=p;};
    void       SetDataPaths(DataPaths p)             { dataPaths_=p; CheckDataPaths(p); };
};

} // close namespace SterileSearch

#endif
