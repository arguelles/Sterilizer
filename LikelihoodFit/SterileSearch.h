#ifndef _H_STERILE_SEARCH_
#define _H_STERILE_SEARCH_


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

class SterileHunter {
  private:
    /// Parameters used in the fit
    double maxFitEnergy=2.e4;
    double minCosth = -1.;
    double maxCosth = 0.2;
    double minAstroEnergy=0.0;
    double maxAstroEnergy=1e10;
    double minAzimuth = 0.;
    double maxAzimuth = 2.*pi<double>();

    // To store best fit point, fit seed, and data challenge
    std::vector<double> existingBest;
    std::vector<double> fitSeed{1.02,0,0,0.05,.0985,1.1,1,0};
    std::vector<double> dataChallenge_nuisance_parameters{1,0,0,0,.1,1,1,0};

    unsigned int number_of_data_challenges = 1;
    // options
    bool use_factorization_technique=false;
    bool use_datachallenge_histogram=false;
    bool writeCompact=false;
    bool readCompact=false;
    bool doDataChallenge=false;
    bool exitAfterLoading=false;
    bool UseBurnsample=true;
    bool use_gsl_minimizer = false;
    bool do_asimov_sensitivity = false;

    int yearsIC86=1;
    // sterile neutrino parameters
    double th24_null = 0;
    double dm41sq_null = 0;
    bool dump_data = false;
    bool dump_mc_data = false;
    bool dump_real_data = false;
    bool dump_fit = false;
    bool save_flux = false;
    bool save_dc_sample = false;

  public:
    /// \brief Constructor
    SterileHunter(){

    };
};

#endif
