#ifndef _H_OVERSIZEWEIGHT_
#define _H_OVERSIZEWEIGHT_

#include "oversizeWeight.h"

//Constructor 
OversizeWeighter::OversizeWeighter(std::string FileName) {
  std::ifstream file(FileName);
  if(!file.good())
    throw std::runtime_error("Oversize file "+FileName+" does not exist!");
  std::string line;
  std::vector<std::string> elements;

  // Get energy bin edges (1st line)
  getline(file,line);
  boost::split(elements, line, boost::is_any_of("\t \n"));    
  for(int i=0; i!=elements.size(); ++i)
      EnergyBinEdges.push_back(stod(elements.at(i)));
  elements.clear();

  // Get zenith bin edges (2nd line)
  getline(file,line);
  boost::split(elements, line, boost::is_any_of("\t "));    
  for(int i=0; i!=elements.size(); ++i)
      ZenithBinEdges.push_back(stod(elements.at(i)));
  elements.clear();

  // Get correction function
  for(int i=0; i!=EnergyBinEdges.size()-1; ++i)
    {
      std::vector<double> EnergyRow;
      getline(file,line);
      boost::split(elements, line, boost::is_any_of("\t "));    
      for(int i=0; i!=elements.size(); ++i)
	{
	  EnergyRow.push_back(stod(elements.at(i)));
	}
      if(EnergyRow.size()==(ZenithBinEdges.size()-1))
	CorrectionFunction.push_back(EnergyRow);
      else
	throw std::runtime_error("Mismatching row length in OversizeWeight!");
    }
}

//Evaluator
double OversizeWeighter::EvaluateOversizeCorrection(double energy, double zenith) const {
  double coszenith=cos(zenith);

  if( ( energy>EnergyBinEdges.at(EnergyBinEdges.size()-1) ) 
      || ( energy<EnergyBinEdges.at(0) ) ) return 0;
  if( ( coszenith>ZenithBinEdges.at(ZenithBinEdges.size()-1) ) 
      || ( coszenith<ZenithBinEdges.at(0) ) ) return 0;

  int ebin=0, zbin=0;
  for(int i=0; i!=EnergyBinEdges.size()-1; ++i)
    if ( (energy > EnergyBinEdges[i]) && (energy < EnergyBinEdges[i+1]) )
      {
	ebin=i;
	break;
      }
  for(int j=0; j!=ZenithBinEdges.size()-1; ++j)
    if ( (coszenith > ZenithBinEdges[j]) && (coszenith < ZenithBinEdges[j+1]) )
      {
	zbin=j;
	break;
      }

  return CorrectionFunction.at(ebin).at(zbin);
}


#endif
