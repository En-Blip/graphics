#pragma once

#include "../../vectors.h"
#include "../../error.h"
#include "../include/readPNG.h"
#include "../include/safeRead.h"
#include <vector>

struct IDHR_t {
    uint2 imgDims; // 4 bytes of width, 4 bytes of height
    char bit_depth, color_type, compression_method, filter_method, interlace_method; 
};

std::vector<color> loadImgFromPng(const std::string& filename, int2& dims);

// helper functions
bool readPNGheader(std::ifstream& file);
bool readIHDR(std::ifstream& file, IDHR_t& settings);
