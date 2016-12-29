#ifndef RUN_SPEC
#define RUN_SPEC

#include <LeptonWeighter/event.h>
#include <LeptonWeighter/particleType.h>
#include <LeptonWeighter/lepton_weighter.h>

struct run{
	const std::string path;
  const std::string filename;
  const particleType pt;
  const double unshadowedFraction ;// aka domefficiency
  LW::SimDetails details;
	run(const std::string& p,const std::string& filename, particleType pt, double unshadowedFraction, LW::SimDetails d):
	path(p),filename(filename),pt(pt),unshadowedFraction(unshadowedFraction),details(d){}
};

#endif
