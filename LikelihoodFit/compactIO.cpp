#include <fstream>
#include <boost/crc.hpp>

#include "compactIO.h"
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
