#ifndef COMPACTIO_H
#define COMPACTIO_H

#include <deque>
#include "Event.h"
#include "runspec.h"

uint32_t getFileChecksum(const std::string& filename);
void splatData(const std::string& filename, const uint32_t progChecksum, const std::deque<Event>& exp, const std::deque<Event>& sim);
void unsplatData(const std::string& filename, const uint32_t expectedChecksum, std::deque<Event>& exp, std::deque<Event>& sim);

#endif //COMPACTIO_H
