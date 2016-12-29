#ifndef OVERCORR_H
#define OVERCORR_H

#include <iostream>
#include <set>

#include <boost/math/constants/constants.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>

#include "tableio.h"
#include "particleType.h"


class OversizeWeighter
{
 public:
  OversizeWeighter(std::string Filename);
  double EvaluateOversizeCorrection(double energy, double zenith) const;
 private:
  std::vector<double> EnergyBinEdges;
  std::vector<double> ZenithBinEdges;
  std::vector<std::vector<double> > CorrectionFunction;
};
  

OversizeWeighter::OversizeWeighter(std::string FileName)
{
  std::ifstream file(FileName);
  std::string line;
  std::vector<std::string> elements;


  // Get energy bin edges (1st line)
  getline(file,line);
  boost::split(elements, line, boost::is_any_of("\t \n"));    
  for(int i=0; i!=elements.size(); ++i)
    {
      EnergyBinEdges.push_back(stod(elements.at(i)));
    }
  elements.clear();

  // Get zenith bin edges (2nd line)
  getline(file,line);
  boost::split(elements, line, boost::is_any_of("\t "));    
  for(int i=0; i!=elements.size(); ++i)
    {
      ZenithBinEdges.push_back(stod(elements.at(i)));
    }
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
	{
	  std::cout<<"Mismatching row length in OversizeWeight!";
	  assert(0);
	}   
    }
  bool verbose=false;
  if(verbose)
    {     
      std::cout<<"Oversize weighter has loaded correction function. "<<std::endl;
      std::cout<<"Energy bins : " ;
      for(int i=0; i!=EnergyBinEdges.size(); ++i)
	std::cout << EnergyBinEdges.at(i)<<", ";
      std::cout<<std::endl<<"Zenith bins : " ;
      for(int i=0; i!=ZenithBinEdges.size(); ++i)
	std::cout << ZenithBinEdges.at(i)<<", ";
      std::cout<<std::endl<<"Contents : " <<std::endl; 
      for(int i=0; i!=EnergyBinEdges.size()-1; ++i)
	{
	  for(int j=0; j!=ZenithBinEdges.size()-1; ++j)
	    std::cout<<CorrectionFunction[i][j]<< " " ;
	  std::cout<<std::endl;
	}
    }
}
 

double  OversizeWeighter::EvaluateOversizeCorrection(double energy, double zenith) const
{
  zenith=cos(zenith);
  //  std::cout<<"energy is " << energy << ", zenith is " << zenith<<std::endl;
  //  std::cout<<"Limits : "<<EnergyBinEdges.at(EnergyBinEdges.size()-1) << ", " << EnergyBinEdges.at(0) << ", " << ZenithBinEdges.at(ZenithBinEdges.size()-1) <<", " << ZenithBinEdges.at(0)<<std::endl;
  
  if( ( energy>EnergyBinEdges.at(EnergyBinEdges.size()-1) ) 
      || ( energy<EnergyBinEdges.at(0) ) ) return 0;
  if( ( zenith>ZenithBinEdges.at(ZenithBinEdges.size()-1) ) 
      || ( zenith<ZenithBinEdges.at(0) ) ) return 0;
  
  int ebin=0, zbin=0;
  for(int i=0; i!=EnergyBinEdges.size()-1; ++i)
    if ( (energy > EnergyBinEdges[i]) && (energy < EnergyBinEdges[i+1]) )
      {
	ebin=i;
	break;
      }
  for(int j=0; j!=ZenithBinEdges.size()-1; ++j)
    if ( (zenith > ZenithBinEdges[j]) && (zenith < ZenithBinEdges[j+1]) )
      {
	zbin=j;
	break;
      }
  
  return CorrectionFunction.at(ebin).at(zbin);
}

#endif //OVERCORR_H
