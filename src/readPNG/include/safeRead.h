#pragma once 
#include "../../error.h"
#include <sstream>
#include <fstream>

int readInt(std::ifstream& file, unsigned short int bytes);
bool safeRead(std::ifstream& file, void* data, uint32_t size);
