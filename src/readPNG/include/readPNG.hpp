#pragma once

#include "../../vectors.h"
#include "../../error.h"
#include "../include/safeRead.h"
#include "datastream.hpp"
#include <cstdint>
#include <vector>

struct IHDR_t {
  uint2 imgDims; // 4 bytes of width, 4 bytes of height
  char bit_depth, color_type, compression_method, filter_method, interlace_method; 
};

struct ZLIB_hdr {
  uint8_t CM;
  uint8_t CINFO;

  uint8_t FLEVEL;
  uint8_t FDICT;
  uint8_t FCHECK;

  uint32_t DICT;
};

std::vector<color> loadImgFromPng(const std::string& filename, int2& dims);

// helper functions
bool readPNGheader(std::ifstream& file);
bool readIHDR(std::ifstream& file, IHDR_t& settings);
bool readZLIB_hdr(datastream& zlib_ds, uint64_t& bit_idx, ZLIB_hdr& header);
bool readUncompressedSet(datastream& zlib_ds, uint64_t& bit_idx, uint8_t*& output_cursor);

/* Filtering algorithms */
uint8_t filter_none(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x);
uint8_t filter_sub(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x);
uint8_t filter_up(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x);
uint8_t filter_avg(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x);
uint8_t filter_paeth(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x);


void read_grayscale(std::vector<color>& color_arr, std::unique_ptr<uint8_t[]>& byte_arr, IHDR_t& settings);
void read_rgb(std::vector<color>& color_arr, std::unique_ptr<uint8_t[]>& byte_arr, IHDR_t& settings);
