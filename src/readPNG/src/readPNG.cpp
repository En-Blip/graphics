#include "../include/readPNG.hpp"
#include "../include/datastream.hpp"
#include "../include/huffmanTree.hpp"
#include "../../error.h"
#include "../include/safeRead.h"
#include <algorithm>
#include <bitset>
#include <cstdint>
#include <fstream>
#include <ios>
#include <ostream>
#include <vector>

/* IMPROVEMENTS
 * instead of getting all the chunks at once, we can use a 32768 byte sliding window
 * (or whatever is defined) and discard datastream values further than the sliding window back
 * This could be done by making the datastream a circular buffer
 * The filter algorithms need a full rework
 *
 */

// see reference in https://www.rfc-editor.org/rfc/rfc1950,
// https://www.rfc-editor.org/rfc/rfc1951#section-2,
// or https://stackoverflow.com/questions/73779035/how-does-zlib-decompression-work-on-a-png-idat-chunk for the summary
std::vector<color> loadImgFromPng(const std::string& filename, int2& dims){
  const std::vector<color> nullvec; // some vector to return upon failure

  // load image and return the first eight bytes
  // open file
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    ERROR("Error opening the file: " + filename);
    return nullvec;
  }

  // read header
  if(!readPNGheader(file)){
    ERROR("PNG header mismatch");
    return nullvec;
  }

  // find and read IHDR chunk to check settings
  IHDR_t settings;
  if(!readIHDR(file, settings)){
    ERROR("No IHDR chunk");
    return nullvec;
  }

  disp(settings.imgDims.x);
  disp(settings.imgDims.y);
  disp(int(settings.bit_depth));
  disp(int(settings.color_type));

  dims.x = settings.imgDims.x;
  dims.y = settings.imgDims.y;

  // read the CRC (and ignore it for now)
  uint32_t CRC = readInt(file, 4);

  // create the datastream object to begin decoding the data
  datastream zlib_ds = datastream(file);
  // create some iterator to keep track of working file idx
  uint64_t bit_idx = 0;

  // read zlib header from input stream
  ZLIB_hdr ds_hdr;
  if(!readZLIB_hdr(zlib_ds, bit_idx, ds_hdr)){
    ERROR("Something wrong with datastream header");
    return nullvec;
  }

  // decoding specific variables
  bool lastBlock = false;
  std::unique_ptr<uint8_t[]> output_ds = std::make_unique<uint8_t[]>(10 * (settings.imgDims.x * 3 + 1) * settings.imgDims.y);
  uint8_t* output_cursor = output_ds.get();

  // read and decode compressed data
  do{
    // read block header from input stream
    uint8_t HDR = zlib_ds.readBits(bit_idx, 3);
    const uint8_t BFINAL = (HDR & 0b100) >> 2;
    const uint8_t BTYPE = ((HDR & 0b1) << 1) | ((HDR & 0b10) >> 1);
#ifdef CHUNKS
    disp(int(BTYPE));
#endif

    assert(BTYPE != 0b11);

    lastBlock = BFINAL;

    // if stored with no compression
    if(BTYPE == 0b00){
      if(!readUncompressedSet(zlib_ds, bit_idx, output_cursor)){
        ERROR("failed to read uncompressed data");
        ASSERT(false);
      }
    }else{
      huffmanTree huffmanCodes;
      huffmanTree distTree;

      uint8_t bitval, extraBits;
      uint32_t bck_dist;
      uint32_t length;
      uint32_t value;

      // if compressed with dynamic Huffman codes
      if(BTYPE == 0b10){
        // read the various lengths of the block values
        uint16_t HLIT = zlib_ds.readBits(bit_idx, 5, true) + 257;
        uint8_t HDIST = zlib_ds.readBits(bit_idx, 5, true) + 1;
        uint8_t HCLEN = zlib_ds.readBits(bit_idx, 4, true) + 4;

        uint8_t code_idx[] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
        int codeCodeLens[19] = { 0 };

        // create the huffman tree to decode huffman codes (huffman code inception)
        uint64_t temp_idx = bit_idx;
        for(uint8_t i = 0; i < HCLEN; i++) 
          codeCodeLens[code_idx[i]] = zlib_ds.readBits(bit_idx, 3, true);

        if(HDIST == 1 && !codeCodeLens[0]) ERROR("Check functionality of special case ");

        huffmanTree codeCodeTree(codeCodeLens, 19);
        // codeCodeTree.printTree();

        int previousCodeLength = -1;
        auto codeLens = std::make_unique<int[]>( HLIT );
        auto distLens = std::make_unique<int[]>( HDIST );

        uint8_t repeatNum;
        for (uint16_t i = 0; i < (HLIT);) {
          value = codeCodeTree.decode(zlib_ds, bit_idx);

          if(value < 16){
            previousCodeLength = value;
            codeLens[i] = previousCodeLength;
            i++;
            continue;
          }else if(value == 16){
            assert(previousCodeLength >= 0);
            repeatNum = 3 + zlib_ds.readBits(bit_idx, 2, true);
            for(uint8_t x = 0; x < repeatNum; x++){
              codeLens[i] = previousCodeLength;
              i++;
            }
            continue;
          }else if(value == 17){
            repeatNum = 3 + zlib_ds.readBits(bit_idx, 3, true);
          }else if(value == 18){ 
            repeatNum = 11 + zlib_ds.readBits(bit_idx, 7, true);
          }
          for(uint8_t x = 0; x < repeatNum; x++){
            codeLens[i] = 0;
            i++;
          }
        }

        for (uint16_t i = 0; i < (HDIST);) {
          value = codeCodeTree.decode(zlib_ds, bit_idx);

          if(value < 16){
            previousCodeLength = value;
            distLens[i] = value;
            i++;
          }else if(value == 16){
            uint8_t repeatNum = 3 + zlib_ds.readBits(bit_idx, 2, true);
            for(uint8_t x = 0; x < repeatNum; x++){
              distLens[i] = previousCodeLength;
              i++;
            }
          }else if(value == 17){
            uint8_t repeatNum = 3 + zlib_ds.readBits(bit_idx, 3, true);
            for(uint8_t x = 0; x < repeatNum; x++){
              distLens[i] = 0;
              i++;
            }
          }else if(value == 18){
            uint8_t repeatNum = 11 + zlib_ds.readBits(bit_idx, 7, true);
            for(uint8_t x = 0; x < repeatNum; x++){
              distLens[i] = 0;
              i++;
            }
          }
        }

#ifdef VERBOSE
        disp("HLIT");
        for(int i = 0; i < HLIT; i++){
          disp(i);
          disp(codeLens.get()[i]);
        }

        disp("HDIST");
        for(int i = 0; i < HDIST; i++){
          disp(i);
          disp(distLens.get()[i]);
        }
#endif

        huffmanCodes.initialize(codeLens.get(), HLIT);
        distTree.initialize(distLens.get(), HDIST);

      }else{
        int codeLens[] = {8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8};

        // int codeLens[288]; // same thing as above but more readable
        // for(uint16_t i = 0; i < 144; i++) codeLens[i] = 8;
        // for(uint16_t i = 144; i < 256; i++) codeLens[i] = 9;
        // for(uint16_t i = 256; i < 280; i++) codeLens[i] = 7;
        // for(uint16_t i = 280; i < 288; i++) codeLens[i] = 8;

        int distLens[] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

        huffmanCodes.initialize(codeLens, sizeof(codeLens) / sizeof(codeLens[0]));
        distTree.initialize(distLens, sizeof(distLens) / sizeof(distLens[0]));
      }

#ifdef VERBOSE          
      bool last_lit = false;
#endif
      while(1){
        value = huffmanCodes.decode(zlib_ds, bit_idx);

        ASSERT(value <= 285); // max value range

        if(value < 256){
          // copy value (literal byte) to output stream
          assert((1 + output_cursor - output_ds.get()) < (settings.imgDims.x * 3 + 1) * settings.imgDims.y + 1);
          *output_cursor = (uint8_t)value;
          output_cursor++;

#ifdef VERBOSE          
          if(!last_lit) std::cout << "\nliteral";
          std::cout << " " << value;
          last_lit = true;
#endif
        }else{
#ifdef VERBOSE          
          last_lit = false;
#endif
          if(value == 256) break; // end of block

          // use value and extra bits to decode true length
          if(value == 285) length = 258;
          else if (value < 265) length = value - 254;
          else{
            uint8_t minLen[] = {3, 4, 5, 6, 7, 8, 9 , 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227};
            extraBits = (value - 261) / 4;

            length = minLen[value-257] + zlib_ds.readBits(bit_idx, extraBits, true);
          }

          // decode distance from input stream
          uint8_t distCode = distTree.decode(zlib_ds, bit_idx);

          // use dist value and extra bits to decode the true distance
          if(distCode >= 4){
            uint32_t minDist[] = {1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577}; // the minimum distance for each code
            extraBits = (distCode - 2) / 2;

            if(extraBits > 8){
              uint8_t temp = zlib_ds.readBits(bit_idx, 8, true);
              extraBits -= 8;
              bck_dist = static_cast<uint32_t>(zlib_ds.readBits(bit_idx, extraBits, true));
              bck_dist <<= 8;
              bck_dist |= temp;
            }else{
              bck_dist = static_cast<uint32_t>(zlib_ds.readBits(bit_idx, extraBits, true));
            }
            bck_dist += minDist[distCode];
          }else{
            bck_dist = distCode + 1;
          }

          // move backwards distance bytes in the output
          // stream, and copy length bytes from this
          // position to the output stream.
#ifdef VERBOSE          
          std::cout << "\nmatch " << length << " " << bck_dist << std::endl; 
#endif
          uint8_t* cpy_cursor = output_cursor - bck_dist;

          assert(cpy_cursor >= output_ds.get());
          assert((length + output_cursor - output_ds.get()) < 10 *(settings.imgDims.x * 3 + 1) * settings.imgDims.y);

          for (uint32_t k = 0; k < length; k++) {
            *output_cursor = *(output_cursor - bck_dist);
            output_cursor++;
          }
        }
      }
    }

  }while(!lastBlock);

  // read ADLER32 checksum

  // unfilter data
  auto colorBuffer = std::vector<color>(settings.imgDims.x * settings.imgDims.y);

  switch(settings.color_type){
    case 0:
      read_grayscale(colorBuffer, output_ds, settings);

      break;
    case 3:
      ERROR("NO PLTE HANDLING IMPLEMENTED");
      assert(false);
      break;
    case 2:
      read_rgb(colorBuffer, output_ds, settings);

      break;
    case 4:
    case 6:
    default:
      ERROR("Invalid color type");
      assert(false);
  }

  // create the color buffer to return
  return colorBuffer;
}

/* algorithms to pack the color buffer depending on the image's color_type */
void read_grayscale(std::vector<color>& color_arr, std::unique_ptr<uint8_t[]>& byte_arr, IHDR_t& settings){
  disp("grayscale");
  uint8_t (*filter[5])(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x) = \
  {filter_none, filter_sub, filter_up, filter_avg, filter_paeth};

  constexpr uint8_t px_size = 1;
  color tmp_col;
  uint8_t* output_cursor = byte_arr.get();
  uint8_t filter_type;
  uint8_t* prev_line = nullptr; // the previous scanline
  uint8_t* cur_line = new uint8_t[settings.imgDims.x*px_size];

  uint8_t filters_used = 0;

  for(int i = 0; i < settings.imgDims.y; i++){
    filter_type = *output_cursor;
    output_cursor++;
    filters_used = filters_used | 1 << filter_type;

    for(int j = 0; j < settings.imgDims.x; j++){
      uint8_t tmp = filter[filter_type](&output_cursor, prev_line, cur_line, px_size, settings, j);
      cur_line[px_size*j] = tmp;
      tmp_col.r = tmp;
      tmp_col.g = tmp;
      tmp_col.b = tmp;
      color_arr.at(i*settings.imgDims.x + j) = tmp_col;
    }

    // move the cur line up and delete prev line
    if(!prev_line) prev_line = new uint8_t[settings.imgDims.x*px_size];
    std::swap(prev_line, cur_line); 
  }

  disp(std::bitset<8>(filters_used));
  delete[] cur_line;
  delete[] prev_line;
}

void read_rgb(std::vector<color>& color_arr, std::unique_ptr<uint8_t[]>& byte_arr, IHDR_t& settings){
  constexpr uint8_t px_size = 3;
  uint8_t* output_cursor = byte_arr.get();
  color tmp_col;
  uint8_t filter_type;
  uint8_t* prev_line = nullptr; // the previous scanline
  uint8_t* cur_line = new uint8_t[settings.imgDims.x*px_size];

  uint8_t (*filter[5])(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x) = \
  {filter_none, filter_sub, filter_up, filter_avg, filter_paeth};

  uint8_t filters_used = 0;

  for(int i = 0; i < settings.imgDims.y; i++){
    filter_type = *output_cursor;
    filters_used = filters_used | 1 << filter_type;
    output_cursor++;

    for(int j = 0; j < settings.imgDims.x; j++){
      tmp_col.r = filter[filter_type](&output_cursor, prev_line, cur_line, px_size, settings, px_size*j);
      cur_line[px_size*j] = tmp_col.r;
      tmp_col.g = filter[filter_type](&output_cursor, prev_line, cur_line, px_size, settings, px_size*j+1);
      cur_line[px_size*j+1] = tmp_col.g;
      tmp_col.b = filter[filter_type](&output_cursor, prev_line, cur_line, px_size, settings, px_size*j+2);
      cur_line[px_size*j+2] = tmp_col.b;
      color_arr.at(i*settings.imgDims.x + j) = tmp_col;
    }

    // move the cur line up and delete prev line
    if(!prev_line) prev_line = new uint8_t[settings.imgDims.x*px_size];
    std::swap(prev_line, cur_line);
  }

  disp(std::bitset<8>(filters_used));
  delete[] cur_line;
  delete[] prev_line;
}

/* Filtering Algorithms */
// for some reason I named the cursor of the decoded 
// until refactoring, im just putting the number of bytes per pixel for color type

/**
 * cursor: the decoded unfilterd bytes
 * prev_line: (Prior(x))
 * cur_line: (Raw(x))
 * pixel_bytesize: bpp
 * settings: for dims
 */
uint8_t filter_none(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x){
  uint8_t return_col = **cursor;
  (*cursor)++;
  return return_col;
}

uint8_t filter_sub(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x){
  // Raw(x) = Sub(x) + Raw(x-bpp)
  uint8_t return_col; 

  if(x < pixel_bytesize) return_col = **cursor;
  else return_col = **cursor + *(cur_line + x - pixel_bytesize);
  (*cursor)++;
  return return_col;
}

uint8_t filter_up(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x){
  // Raw(x) = Up(x) + Prior(x)
  uint8_t return_col; 

  uint8_t prior_x = (prev_line) ? *(prev_line + x) : 0;
  return_col = **cursor + prior_x;
  (*cursor)++;
  return return_col;
}

uint8_t filter_avg(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x){
  // Raw(x) = Average(x) + floor((Raw(x-bpp)+Prior(x))/2)
  uint8_t return_col;

  uint8_t raw_back = (x >= pixel_bytesize) ? *(cur_line + x - pixel_bytesize) : 0;
  uint8_t prior_x = (prev_line) ? *(prev_line + x) : 0;
  uint16_t avg = ((uint16_t)raw_back + (uint16_t)prior_x)/2; // casting to uin16_t to prevent overflow
  return_col = **cursor + (uint8_t)avg;
  (*cursor)++;
  return return_col;
}

static uint8_t PaethPredictor(uint8_t a, uint8_t b, uint8_t c){
  int p = (int)a + (int)b - (int)c;
  int pa = abs(p - a);
  int pb = abs(p - b);
  int pc = abs(p - c);

  if(pa <= pb && pa <= pc) return a;
  else if(pb <= pc) return b;
  else return c;
}

uint8_t filter_paeth(uint8_t** cursor, uint8_t* prev_line, uint8_t* cur_line, uint8_t pixel_bytesize, IHDR_t& settings, int x){
  // Raw(x) = Paeth(x) + PaethPredictor(Raw(x-bpp), Prior(x), Prior(x-bpp))
  uint8_t return_col;

  uint8_t raw_back = (x >= pixel_bytesize) ? *(cur_line + x - pixel_bytesize) : 0;
  uint8_t prior_x = (prev_line) ? *(prev_line + x) : 0;
  uint8_t prior_back = ((x >= pixel_bytesize) && (prev_line)) ? *(prev_line + x - pixel_bytesize) : 0;

  return_col = **cursor + PaethPredictor(raw_back, prior_x, prior_back);
  (*cursor)++;
  return return_col;
}


bool readPNGheader(std::ifstream& file){
  std::unique_ptr<char[]> header = std::make_unique<char[]>(8);
  char PNGheader[9] = { -119, 80, 78, 71, 13, 10, 26, 10, 0 }; // -119 is 137 as a char
  safeRead(file, header.get(), 8);
  for (int i = 0; i < 8; ++i) {
    if (header[i] != PNGheader[i]) {
      std::cout << header << " != " << PNGheader << std::endl;
      return false;
    }
  }
  return true;
}

bool readIHDR(std::ifstream& file, IHDR_t& settings){
  uint32_t chunkSize;
  safeRead(file, &chunkSize, 4);

  // read the chunk type, since its the first chunk it should be IHDR
  char chunkType[4];
  safeRead(file, chunkType, 4);
  if (chunkType[0] != 'I' || chunkType[1] != 'H' || chunkType[2] != 'D' || chunkType[3] != 'R') {
    return false;
  }

  settings.imgDims.x = readInt(file, 4);
  settings.imgDims.y = readInt(file, 4); 

  // if this is IHDR, will contain, in bytes
  // width(4), height(4), bit depth(1), color type(1), compression method(1)
  // filter method(1), interlace method(1)
  safeRead(file, &settings.bit_depth, 5); // careful: cpp compiler inserts more bited in struct
                                          // for alignment -> dont use sizeof

                                          // make sure the settings work for what I have 
                                          // ASSERT(settings.imgDims.x * settings.imgDims.y == 3024 * 1964);
                                          // ASSERT(settings.color_type == 2);
                                          // ASSERT(settings.bit_depth == 8 || settings.bit_depth == 16);
  disp(int(settings.color_type));
  ASSERT(settings.bit_depth == 8);
  ASSERT(settings.compression_method == 0);
  ASSERT(settings.filter_method == 0);
  ASSERT(settings.interlace_method == 0);

  return true;
}

bool readZLIB_hdr(datastream& zlib_ds, uint64_t& bit_idx, ZLIB_hdr& header) {

  uint8_t CMF = zlib_ds.readBits(bit_idx, 8, true);
  header.CM = CMF & 0x0F; // compression method
  header.CINFO = (CMF & 0xF0) >> 4; // compression info 

  ASSERT(header.CM == 8); // CM should be 8 for PNG files
  ASSERT(header.CINFO < 8); // values above 7 not allows in this version of the specification

  // uint32_t wind_sz = 1u << (CINFO + 8); // or is it pow(2, cinfo) + 8; ??
  // window size doesnt really matter for us because I
  // want to limit file reads anyways, but may be a 
  // good redundancy check in the future

  uint8_t FLG =  zlib_ds.readBits(bit_idx, 8, true);
  header.FLEVEL = (FLG & 0xc0) >> 6;
  header.FDICT = (FLG & 0x20) >> 5;
  header.FCHECK = FLG & 0x1F;

  assert(((CMF << 8) + FLG) % 31 == 0); // check for validity of CMF and FLG;

  // read dict if required
  if(header.FDICT){
    header.DICT = zlib_ds.readInt_ds(bit_idx);
  }

  return true;
}

bool readUncompressedSet(datastream& zlib_ds, uint64_t& bit_idx, uint8_t*& output_cursor){
  // if stored with no compression
  // skip any remaining bits in current partially processed byte
  bit_idx += (bit_idx % 8) ? (8 - bit_idx % 8) : 0;

  // read LEN and NLEN
  uint16_t LEN = zlib_ds.readInt_ds(bit_idx, 2);
  uint16_t NLEN = zlib_ds.readInt_ds(bit_idx, 2);

  ASSERT((uint16_t)(~LEN) == NLEN);

  // copy LEN bytes of data to output
  while(LEN > 0){
    uint32_t bytes_read = std::min(LEN, static_cast<uint16_t>(zlib_ds.wind_sz));
    zlib_ds.readBytes(output_cursor, bit_idx, bytes_read); // LEN will never be more 

    LEN -= bytes_read;
    output_cursor += bytes_read;
  }

  return true;
}
