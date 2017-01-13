#ifndef OVERCORR_H
#define OVERCORR_H

// This class exists to apply an event weight read from a binned table to data or MC.
//  We use it to correct oversized MC to non-oversized distributions,

#include <iostream>
#include <set>

#include <boost/math/constants/constants.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <string>
#include <vector>

class OversizeWeighter
{
 public:
  // Constructor
  OversizeWeighter(std::string Filename);
  OversizeWeighter() {};
  // Evaluator
  double EvaluateOversizeCorrection(double energy, double zenith) const;
 private:
  std::vector<double> EnergyBinEdges;
  std::vector<double> ZenithBinEdges;
  std::vector<std::vector<double> > CorrectionFunction;
};

#endif //OVERCORR_H
