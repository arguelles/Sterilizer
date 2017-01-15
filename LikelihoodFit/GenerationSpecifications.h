#ifndef _H_GENERATION_SPECIFICATIONS_
#define _H_GENERATION_SPECIFICATIONS_

#include <LeptonWeighter/lepton_weighter.h>
#include <LeptonWeighter/event.h>
#include <LeptonWeighter/particleType.h>
#include "runspec.h"
using boost::math::constants::pi;

std::map<std::string,run> GetSimInfo(std::string mc_dataPath)
{
   std::map<std::string,run> simInfo = {
                                {"nufsgen_mie_0_9",run(mc_dataPath,"noover_mie_0.9.h5",particleType::NuMu,0.9,
                                        LW::SimDetails {1,615360000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_mie_0_95",run(mc_dataPath,"noover_mie_0.95.h5",particleType::NuMu,0.95,
                                        LW::SimDetails {1,292140000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_mie_0_99",run(mc_dataPath,"noover_mie_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,734400000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_mie_1_089",run(mc_dataPath,"noover_mie_1.089.h5",particleType::NuMu,1.089,
                                        LW::SimDetails {1,584640000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_mie_1_1979",run(mc_dataPath,"noover_mie_1.1979.h5",particleType::NuMu,1.1979,
                                        LW::SimDetails {1,440400000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_mie_0_99_sel",run(mc_dataPath,"noover_mie_0.99_sel.h5",particleType::NuMu,0.99,
							LW::SimDetails {1,57300000, //total events
							2011,//year
							800, //injection radius
							0,2*pi<double>(), //azimuth range
							80.0*pi<double>()/180.0,pi<double>(), //zenith range
							2e2,1e6, //energy range
							2} // spectral index
				)},
                                {"nufsgen_lea_0_9",run(mc_dataPath,"noover_lea_0.9.h5",particleType::NuMu,0.9,
                                        LW::SimDetails {1,39416000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_lea_0_99",run(mc_dataPath,"noover_lea_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,39568000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_lea_1_1979",run(mc_dataPath,"noover_lea_1.1979.h5",particleType::NuMu,1.1979,
                                        LW::SimDetails {1,39616000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_noholeice_mie_0_99",run(mc_dataPath,"noholeice_mie_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,78856000, //total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_ice1_mie_0_99",run(mc_dataPath,"ice1_mie_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,79536000,//total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )},
                                {"nufsgen_ice2_mie_0_99",run(mc_dataPath,"ice2_mie_0.99.h5",particleType::NuMu,0.99,
                                        LW::SimDetails {1,79584000,//total events
                                                        2011,//year
                                                        800, //injection radius
                                                        0,2*pi<double>(), //azimuth range
                                                        80.0*pi<double>()/180.0,pi<double>(), //zenith range
                                                        2e2,1e6, //energy range
                                                        2} // spectral index
                                )}
                              };

  return simInfo;
}
#endif
