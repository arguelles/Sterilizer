#ifndef _H_STERILE_SEARCH_
#define _H_STERILE_SEARCH_

struct fitResult {
  std::vector<double> params;
  double likelihood;
  unsigned int nEval, nGrad;
  bool succeeded;
};

struct SterileNeutrinoParameters {
  double th14 = 0;
  double th24 = 0;
  double th34 = 0;
  double del14 = 0;
  double del24 = 0;
  double dm41sq = 0;
  SterileNeutrinoParameters(double th14,double th24, double th34, double del14, double del24, double dm41sq):
    th14(th14),th24(th24),th34(th34),del14(del14),del24(del24)
};

using namespace phys_tools::histograms;
using namespace likelihood;
using HistType = histogram<3,entryStoringBin<std::reference_wrapper<const Event>>>;

auto binner = [](HistType& h, const Event& e){
  h.add(e.energy,cos(e.zenith),e.year,amount(std::cref(e)));
};

class SterileHunter {
  private:
    std::vector<unsigned int> years;
    /// Parameters used in the fit
    double maxFitEnergy_ = 2.e4;
    double minCosth_ = -1.;
    double maxCosth_ = 0.2;
    double minAstroEnergy_ =0.0;
    double maxAstroEnergy_ =1e10;
    double minAzimuth_ = 0.;
    double maxAzimuth_ = 2.*pi<double>();

    // To store best fit point, fit seed, and data challenge
    std::vector<double> existingBest_;
    std::vector<double> fitSeed_{1.02,0,0,0.05,.0985,1.1,1,0};
    std::vector<double> dataChallenge_nuisance_parameters_{1,0,0,0,.1,1,1,0};

    unsigned int number_of_data_challenges_ = 1;
    // options
    bool use_factorization_technique_=false;
    bool use_datachallenge_histogram_=false;
    bool writeCompact_=false;
    bool readCompact_=false;
    bool doDataChallenge_=false;
    bool exitAfterLoading_=false;
    bool UseBurnsample_=true;
    bool use_gsl_minimizer_ = false;
    bool do_asimov_sensitivity_ = false;

    int yearsIC86_=1;
    // sterile neutrino parameters
    double th24_null_ = 0;
    double dm41sq_null_ = 0;
    bool dump_data_ = false;
    bool dump_mc_data_ = false;
    bool dump_real_data_ = false;
    bool dump_fit_ = false;
    bool save_flux_ = false;
    bool save_dc_sample_ = false;

    std::vector<double> livetime_;

    // to store events
    std::deque<Event> mainSimulation;
    std::deque<Event> alternativeSimulation;
    std::deque<Event> sample;

    // random number generator
    unsigned int rng_seed;
    std::mt19937 rng;

    // histograms
    HistType data_hist;
    HistType sim_hist;

    // weighter object
    DiffuseFitWeighterMaker DFWM;
    std::shared_ptr<LW::Flux> flux_kaon,flux_pion,flux_prompt;
    std::shared_ptr<LW::CrossSectionFromSpline> xsw;
    LW::mcgenWeighter mcw;
    LW::LeptonWeighter PionFluxWeighter;
    LW::LeptonWeighter KaonFluxWeighter;
    LW::LeptonWeighter PromptFluxWeighter;

    // DOM efficiency splines
    std::vectior<std::unique_ptr<Splinetable>> domEffConv;
    std::vectior<std::unique_ptr<Splinetable>> domEffPrompt;
  public:
    /// \brief Constructor
    SterileHunter(){
      if(readCompact_){
        LoadCompactData(data_filepath);
        LoadCompactMC(mc_filepath);
      } else {
        LoadData(data_filepath);
        LoadMC(mc_filepath);
      }
      LoadFluxes(flux_path,snp)
    }
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
    void ConstructFluxWeighter(std::string squids_files_path,std::string splines_path,SterileNeutrinoParameters snp);
  public:
    marray<double,3> GetDataDistribution();
    marray<double,3> GetExpectation(SterileNeutrinoParameters snp, std::vector<double> nuisance);
    marray<double,3> GetRealization(SterileNeutrinoParameters snp, std::vector<double> nuisance);
    double llhFull(SterileNeutrinoParameters snp, std::vector<double> nuisance){}
    fitResult llh(SterileNeutrinoParameters snp) {}
    // set functions
    void SetRandomNumberGeneratorSeed(unsigned int seed) { }
    unsigned int GetRandomNumberGeneratorSeed(unsigned int seed) { return(rng_seed); }
    void SetUseBurnSample(bool ubs) { UseBurnsample = ubs; }
    bool GetUseBurnSample() { return(UseBurnsample); }
};

#endif
