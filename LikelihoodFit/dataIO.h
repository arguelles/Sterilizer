#ifndef DATAIO_H
#define DATAIO_H

#include <deque>

#include <Event.h>
#include "runspec.h"

std::deque<Event> loadExperimentalData(const std::string& dataPath, bool UseBurnsample);
void loadSimulatedData(std::deque<Event>& buffer, const std::string& dataPath, const std::map<unsigned int,double>& livetime, const std::map<std::string,run>& simInfo,
                       std::vector<std::string> simSetsToLoad, bool loadTargeted);

uint32_t getFileChecksum(const std::string& filename);
void splatData(const std::string& filename, const uint32_t progChecksum, const std::deque<Event>& exp, const std::deque<Event>& sim);
void unsplatData(const std::string& filename, const uint32_t expectedChecksum, std::deque<Event>& exp, std::deque<Event>& sim);

#endif //DATAIO_H
