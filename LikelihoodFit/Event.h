#ifndef EVENT_H
#define EVENT_H

#include <iostream>
#include <set>

#include <boost/math/constants/constants.hpp>

#include "tableio.h"
#include "particleType.h"


double differenceAngle(double zenith1, double azimuth1, double zenith2, double azimuth2);
herr_t collectTableNames(hid_t group_id, const char * member_name, void* operator_data);

struct vector{
	double x,y,z;
	vector():x(0),y(0),z(0){}
	vector(double x, double y, double z):x(x),y(y),z(z){}
	vector(double zenith, double azimuth){
		using boost::math::constants::pi;
		const double theta = pi<double>()-zenith;
		const double phi = azimuth-pi<double>();
		const double rho = std::sin(theta);
		x = rho*std::cos(phi);
		y = rho*std::sin(phi);
		z = std::cos(theta);
	}
	
	double operator*(const vector& rhs) const{
		return x*rhs.x + y*rhs.y + z*rhs.z;
	}
	vector& operator*=(double a){
		x*=a;
		y*=a;
		z*=a;
		return *this;
	}
	vector operator*(double a) const{
		return vector(*this)*=a;
	}
	vector& operator-=(const vector& rhs){
		x-=rhs.x;
		y-=rhs.y;
		z-=rhs.z;
		return *this;
	}
	vector operator-(const vector& rhs) const{
		return vector(*this)-=rhs;
	}
	double magnitude() const{
		return sqrt(*this * *this);
	}
};

enum class Level {generation,level_1,level_2,neutrino};

class Event{
public:
  // event ids
  unsigned int runid;
  unsigned int eventid;
  // new mc variables
  double leptonEnergyFraction;
  double injectedEnergy;
  double totalColumnDepth;
  double inelasticityProbability;
  double intX,intY;

  // true muon quantities
  double injectedMuonEnergy;
  double injectedMuonZenith;
  double injectedMuonAzimuth;

  // reconstructed quantities
  double zenith;
  double azimuth;
  double energy;

  // cached quantities
  double cachedLivetime;
  double cachedNugenWeight;
  double cachedConvPionWeight;
  double cachedConvKaonWeight;
  double cachedPromptWeight;
  double cachedAstroWeight;
  double cachedWeight;

  // dom eff structure
  struct domEffValues{
		double baseRate;

		domEffValues():baseRate(std::numeric_limits<double>::quiet_NaN()){}
	};
	domEffValues cachedConvDOMEff;
	domEffValues cachedPromptDOMEff;
	domEffValues cachedAstroDOMEff;

  // not pass analysis cuts, i.e. was cutted
  bool cutL3;
  int paraboloidStatus;

  particleType primaryType;
  unsigned int year;
  double livetime;

  Event():
    runid(std::numeric_limits<unsigned int>::quiet_NaN()),
    eventid(std::numeric_limits<unsigned int>::quiet_NaN()),
    energy(std::numeric_limits<float>::quiet_NaN()),
    zenith(std::numeric_limits<float>::quiet_NaN()),
    azimuth(std::numeric_limits<float>::quiet_NaN()),
    injectedEnergy(std::numeric_limits<float>::quiet_NaN()),
    leptonEnergyFraction(std::numeric_limits<float>::quiet_NaN()),
    totalColumnDepth(std::numeric_limits<float>::quiet_NaN()),
    inelasticityProbability(std::numeric_limits<float>::quiet_NaN()),
    intX(std::numeric_limits<float>::quiet_NaN()),
    intY(std::numeric_limits<float>::quiet_NaN()),

    livetime(std::numeric_limits<double>::quiet_NaN()),
    primaryType(particleType::unknown),
    cutL3(false),paraboloidStatus(0),
    year(0)
  {}

  bool check(bool check_, Level level) const{
    std::cout << cutL3 << " " << paraboloidStatus << std::endl;
    if (!check_)
      return true;
    else if ( level == Level::generation )
      return(!std::isnan(injectedEnergy) && !std::isnan(injectedMuonAzimuth) && !std::isnan(injectedMuonZenith)
         && !std::isnan(totalColumnDepth) && !std::isnan(inelasticityProbability)
         && !std::isnan(intX) && !std::isnan(intY)
         && primaryType!=particleType::unknown);
    else if ( level == Level::neutrino ){
      return(!std::isnan(injectedEnergy) &&!std::isnan(injectedMuonAzimuth) && !std::isnan(injectedMuonZenith)
             && !std::isnan(totalColumnDepth) && !std::isnan(inelasticityProbability)
             && !std::isnan(intX) && !std::isnan(intY)
             && !std::isnan(zenith) && !std::isnan(energy) && !std::isnan(azimuth)
             && !cutL3 && paraboloidStatus==0
             && primaryType!=particleType::unknown);
    } else {
      return false;
    }
  }
};

namespace{
  template<typename T>
  using simpleValue = TableRow<field<T,CTS("value")>>;
}

template<typename CallbackType>
void readFile(const std::string& filePath, CallbackType action){
  H5File h5file(filePath);
  if(!h5file)
    throw std::runtime_error("Unable to open "+filePath);
  std::set<std::string> tables;
  H5Giterate(h5file,"/",NULL,&collectTableNames,&tables);
  if(tables.empty())
    throw std::runtime_error(filePath+" contains no tables");
  std::map<RecordID,Event> intermediateData;

  // DO NOT change anything above this line

  if(tables.count("MuEx")){
    using particle = TableRow<field<double,CTS("time")>,
                              field<double,CTS("zenith")>,
                              field<double,CTS("azimuth")>,
                              field<double,CTS("energy")>,
                              field<int,CTS("type")>>;

    readTable<particle>(h5file, "MuEx", intermediateData,
          [](const particle& p, Event& e){
            e.runid = p.id.run;
            e.eventid = p.id.event;
            e.zenith=p.get<CTS("zenith")>();
            e.azimuth=p.get<CTS("azimuth")>();
            e.energy=p.get<CTS("energy")>();
          });
  }

  if(tables.count("GenerationSpec")){
    using nmcprop = TableRow<field<double,CTS("injectedEnergy")>,
                             field<double,CTS("muonEnergyFraction")>,
                             field<double,CTS("finalStateProbability")>,
                             field<double,CTS("finalStateX")>,
                             field<double,CTS("finalStateY")>,
                             field<double,CTS("totalColumnDepth")>>;

    readTable<nmcprop>(h5file,"GenerationSpec",intermediateData,
                 [](const nmcprop& p, Event& e){
                    e.injectedEnergy = p.get<CTS("injectedEnergy")>();
                    e.leptonEnergyFraction = p.get<CTS("muonEnergyFraction")>();
                    e.totalColumnDepth = p.get<CTS("totalColumnDepth")>();
                    e.inelasticityProbability = p.get<CTS("finalStateProbability")>();
                    e.intX = p.get<CTS("finalStateX")>();
                    e.intY = p.get<CTS("finalStateY")>();
                 });
  }

  { // functions in this unnamed focus are auxiliary and will be remove eventually
    using iniceprop = TableRow<field<double,CTS("zenith")>,
                               field<double,CTS("azimuth")>,
                               field<double,CTS("energy")>,
                               field<int,CTS("pdg_encoding")>,
                               field<int,CTS("type")>>;

    if(tables.count("MCMaxInIceTrack"))
      readTable<iniceprop>(h5file,"MCMaxInIceTrack",intermediateData,
                   [](const iniceprop& p, Event& e){
                    e.injectedMuonEnergy=p.get<CTS("energy")>();
                    e.injectedMuonZenith=p.get<CTS("zenith")>();
                    e.injectedMuonAzimuth=p.get<CTS("azimuth")>();
                    e.primaryType=static_cast<particleType>(p.get<CTS("type")>());
                    });

    if(tables.count("CutL3")){
      using cutFlag = TableRow<field<unsigned char,CTS("value")>>;
      readTable<cutFlag>(h5file, "CutL3", intermediateData,
      [](const cutFlag& c, Event& e){ e.cutL3=c.get<CTS("value")>(); });
    }

    if(tables.count("TrackFitParaboloidFitParams")){
      using paraboloidParams = TableRow<field<int,CTS("status")>>;
      readTable<paraboloidParams>(h5file, "TrackFitParaboloidFitParams", intermediateData,
      [](const paraboloidParams& c, Event& e){ e.paraboloidStatus = c.get<CTS("status")>(); });
    }
  }

  // The following lines are here to temporarly accomodate the new file format before the new LW
  // files are ready. Please if you do not like this keep hitting Carlos to make LI/LW work
  // this is ugly. Love. CA.

  {
    if(tables.count("WeightingMuon")){
      using weightingmuon= TableRow<field<double,CTS("zenith")>,
                                    field<double,CTS("azimuth")>,
                                    field<double,CTS("energy")>,
                                    field<int,CTS("type")>>;

      readTable<weightingmuon>(h5file,"WeightingMuon",intermediateData,
                   [](const weightingmuon& p, Event& e){
                    e.injectedMuonEnergy=p.get<CTS("energy")>();
                    e.injectedMuonZenith=p.get<CTS("zenith")>();
                    e.injectedMuonAzimuth=p.get<CTS("azimuth")>();
                    e.primaryType=static_cast<particleType>(p.get<CTS("type")>());
                    });

    }
  }

  // DO NOT modify the following line.

  for(std::map<RecordID,Event>::value_type& item : intermediateData)
    action(item.first,item.second);
}

#endif //EVENT_H
