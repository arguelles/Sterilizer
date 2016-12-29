#include <fstream>
#include <boost/crc.hpp>

#include "dataIO.h"
//#include "targetedSim.h"
#include "analysisWeighting.h"

namespace{
	struct energyRecord{
		double energy;
		RecordID id;
	};
	
	bool operator<(const energyRecord& e1, const energyRecord& e2){
		return(e1.energy<e2.energy);
	}
}

std::deque<Event> loadExperimentalData(const std::string& dataPath, bool UseBurnsample){
	std::deque<Event> storage;
	auto dataAction = [&](RecordID id, Event& e, int dataYear){
		if(e.check(false,Level::neutrino)){
			e.year=dataYear;
			storage.push_back(e);
		}
	};
	
	//auto ic79Action=[&](RecordID id, Event& e){ dataAction(id,e,2010); };
	//readFile(dataPath+"burnsample_ic79.h5",ic79Action);
	//readFile(dataPath+"IC79.h5",ic79Action);
	
	auto ic86Action=[&](RecordID id, Event& e){ dataAction(id,e,2011); };
	//readFile(dataPath+"burnsample_ic86_ereco.h5",ic86Action);
  if (UseBurnsample)
    readFile(dataPath+"burnsample_ic86.h5",ic86Action);
  else
    readFile(dataPath+"IC86.h5",ic86Action);
	return(storage);
}

namespace{
	struct domEffSetter{
		simpleEffRate<Event> convDOMEffRate;
		//simpleEffRate<Event> promptDOMEffRate;
		//simpleEffRate<Event> astroDOMEffRate;
		
		domEffSetter(double simulatedDOMEfficiency):
		convDOMEffRate(domEffConv2011.get(),simulatedDOMEfficiency,&Event::cachedConvDOMEff)//,
		//promptDOMEffRate(domEffPrompt2011.get(),simulatedDOMEfficiency,&Event::cachedPromptDOMEff),
		//astroDOMEffRate(domEffAstro2011.get(),simulatedDOMEfficiency,&Event::cachedAstroDOMEff)
		{}
		
		void setCache(Event& e) const{
			convDOMEffRate.setCache(e);
		//	promptDOMEffRate.setCache(e);
		//	astroDOMEffRate.setCache(e);
		}
	};
}

void loadSimulatedData(std::deque<Event>& buffer, const std::string& dataPath, const std::map<int,double>& livetime, const std::map<std::string,run>& simInfo,
                                    std::vector<std::string> simSetsToLoad, bool loadTargeted){
	auto simAction=[&](RecordID id, Event& e, int simYear, const domEffSetter& domEff){
		if(e.check(true,Level::neutrino) && e.energy>1){
			e.year=simYear;
			e.cachedLivetime=livetime.find(simYear)->second;
			e.cachedConvPionWeight=0;
			e.cachedConvKaonWeight=0;
			e.cachedPromptWeight=0;
			e.cachedAstroWeight=0;
			domEff.setCache(e);
			if(e.primaryType==particleType::NuTau || e.primaryType==particleType::NuTauBar){
				assert(e.cachedConvPionWeight==0.0);
				assert(e.cachedConvKaonWeight==0.0);
				assert(e.cachedPromptWeight==0.0);
			}

			buffer.push_back(e);
		}
	};
	
	for(auto simSet : simSetsToLoad){
		const auto& setInfo=simInfo.find(simSet)->second;
		int simYear=setInfo.details.year;
		domEffSetter domEff(setInfo.unshadowedFraction);
		auto callback=[&,simYear](RecordID id, Event& e){ simAction(id,e,simYear,domEff); };
    auto path=dataPath+setInfo.filename;
    readFile(path,callback);
	}
  /*
	if(loadTargeted){
		std::map<std::string,unsigned int> sets{{"Delta/IC79/eff0.9900",2010},{"Echo/IC86/eff0.9900",2011}};
		const double unshadowedFraction=0.99;
		domEffSetter domEff(unshadowedFraction);
		for(auto set : sets){
			int simYear=set.second;
			auto callback=[&,simYear](RecordID id, Event& e){ simAction(id,e,simYear,domEff); };
			for(unsigned int i=0; i<1100; i+=100)
				readFile(dataPath+set.first+"/"+boost::lexical_cast<std::string>(i)+"-"+boost::lexical_cast<std::string>(i+99)+".h5",callback);
		}
	}
  */
}

uint32_t getFileChecksum(const std::string& filename){
	const std::streamsize bufferSize=1u<<16; //1MB
	
	std::ifstream file(filename, std::ios_base::binary);
	if(!file)
		throw std::runtime_error("Unable to open "+filename+" for reading");
	
	boost::crc_32_type fileCRC;
	char buffer[bufferSize];
	do{
		file.read(buffer, bufferSize);
		fileCRC.process_bytes(buffer, file.gcount());
	} while(file);
	
	return(fileCRC.checksum());
}

void splatData(const std::string& filename, const uint32_t progChecksum, const std::deque<Event>& exp, const std::deque<Event>& sim){
	std::ofstream datafile(filename);
	if(!datafile)
		throw std::runtime_error("Unable to open "+filename+" for writing");
	size_t size;
	boost::crc_32_type fileCRC;
	
	datafile.write((char*)&progChecksum,sizeof(progChecksum));
	fileCRC.process_bytes((char*)&progChecksum,sizeof(progChecksum));
	
	size=exp.size();
	datafile.write((char*)&size,sizeof(size));
	fileCRC.process_bytes((char*)&size,sizeof(size));
	for(const Event& e : exp){
		datafile.write((char*)&e,sizeof(e));
		fileCRC.process_bytes((char*)&e,sizeof(e));
	}
	
	size=sim.size();
	datafile.write((char*)&size,sizeof(size));
	fileCRC.process_bytes((char*)&size,sizeof(size));
	for(const Event& e : sim){
		datafile.write((char*)&e,sizeof(e));
		fileCRC.process_bytes((char*)&e,sizeof(e));
	}
	
	uint32_t checksum=fileCRC.checksum();
	datafile.write((char*)&checksum,sizeof(checksum));
}

void unsplatData(const std::string& filename, const uint32_t expectedChecksum, std::deque<Event>& exp, std::deque<Event>& sim){
	std::ifstream datafile(filename);
	if(!datafile)
		throw std::runtime_error("Unable to open "+filename+" for reading");
	boost::crc_32_type fileCRC;
	
	//read checksum; check it
	uint32_t storedChecksum;
	datafile.read((char*)&storedChecksum,sizeof(storedChecksum));
	fileCRC.process_bytes((char*)&storedChecksum,sizeof(storedChecksum));
	if(storedChecksum!=expectedChecksum){
		std::ostringstream ss;
		ss << "Program checksum stored in " << filename << ", " << std::hex << storedChecksum << ", does not match expected (current) checksum, " << expectedChecksum;
		throw std::runtime_error(ss.str());
	}
	
	size_t size;
	const size_t maxEvents=5e7; //as a vague sort-of safety check assume that there will never be more than 50 million events
	datafile.read((char*)&size,sizeof(size));
	fileCRC.process_bytes((char*)&size,sizeof(size));
	if(size>maxEvents)
		throw std::runtime_error(filename+" claims to contain "+boost::lexical_cast<std::string>(size)
								 +" experimental events, which is larger than the safety limit of "
								 +boost::lexical_cast<std::string>(maxEvents));
	exp.resize(size);
	for(const Event& e : exp){
		datafile.read((char*)&e,sizeof(e));
		fileCRC.process_bytes((char*)&e,sizeof(e));
	}
	
	datafile.read((char*)&size,sizeof(size));
	fileCRC.process_bytes((char*)&size,sizeof(size));
	if(size>maxEvents)
		throw std::runtime_error(filename+" claims to contain "+boost::lexical_cast<std::string>(size)
								 +" simulated events, which is larger than the safety limit of "
								 +boost::lexical_cast<std::string>(maxEvents));
	sim.resize(size);
	for(const Event& e : sim){
		datafile.read((char*)&e,sizeof(e));
		fileCRC.process_bytes((char*)&e,sizeof(e));
	}
	
	datafile.read((char*)&storedChecksum,sizeof(storedChecksum));
	if(storedChecksum!=fileCRC.checksum()){
		std::ostringstream ss;
		ss << filename << " appears to be corrupted: stored checksum, " << std::hex << storedChecksum << ", does not match recomputed checksum, " << fileCRC.checksum();
		throw std::runtime_error(ss.str());
	}
}
