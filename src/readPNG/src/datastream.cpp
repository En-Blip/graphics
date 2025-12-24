#include "../include/datastream.hpp"
#include "../include/safeRead.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>

// TODO: split off a load function from operator[] then use that in readBytes

// DOESNT REVERSE BIT ORDER
uint8_t datastream::operator[](int idx){
function_start: // label to replace recursive call when skipping ancillary chunk
                // not sure if there is a better way to do this other than a big loop

                // convert the index into a multiple of the size of the array
  idx %= wind_sz * 3;

  // find which half the access is requesting    
  int8_t working_half = idx / wind_sz;

  // return an error if we are accessing beyond end_idx
  uint8_t end_working_half = end_idx / wind_sz;
  if(end_idx >= 0 && ((working_half == end_working_half && idx > end_idx) || working_half == (end_working_half + 1)%3)){
    ERROR("Access out of datastream bounds");
    // ASSERT(false);
  }

  // load the appropriate half if needed
  if(loaded_half < working_half || (loaded_half == 2 && working_half == 0)){ // im gonna ignore
                                                                             // this magic number

    uint32_t bytes_read = 0;
    uint32_t bytes_to_read;
    unsigned char* load_point = data + wind_sz * working_half;

    while(bytes_read < wind_sz){
      // read the remaining bytes in chunk
      if(bytes_left_in_chunk){
        bytes_to_read = std::min(bytes_left_in_chunk, wind_sz - bytes_read);
        safeRead(file, load_point + bytes_read, bytes_to_read);
        bytes_read += bytes_to_read;
        bytes_left_in_chunk -= bytes_to_read;

        // read the chunk CRC (and ignore it)
        if(!bytes_left_in_chunk)
          uint32_t CRC = readInt(file, 4);
      }

      if(bytes_read < wind_sz){ 
        uint32_t chunkSize = readInt(file, 4);

        // read the chunk type
        char chunkType[5] = {0};
        safeRead(file, chunkType, 4);

#ifdef CHUNKS
        disp(chunkType);
        disp(chunkSize);
#endif

        // confirm that it is an IDAT chunk
        if(!strcmp(chunkType, "IEND")){
          end_idx = wind_sz * working_half + bytes_read;
          if(working_half == end_working_half && idx > end_idx){
            ERROR("Access out of datastream bounds");
            ASSERT(working_half != end_working_half && idx <= end_idx);
          }
          break;
        }else if(!strcmp(chunkType, "IDAT")){
          // read the data and return it
          bytes_to_read = std::min(chunkSize, wind_sz - bytes_read);
          safeRead(file, load_point + bytes_read, bytes_to_read);
          bytes_read += bytes_to_read;
          bytes_left_in_chunk = chunkSize - bytes_to_read;

          // read the chunk CRC (and ignore it)
          if(!bytes_left_in_chunk)
            uint32_t CRC = readInt(file, 4);

        }else if(chunkType[0] & 0x20){ // ancillary chunk
          std::string type(chunkType);
          ERROR("Skipping chunk type: " + type);
          file.seekg(chunkSize + 4, std::ios::cur);
          goto function_start; // tail recursion
        }else{
          std::string type(chunkType);
          ERROR("Unknown chunk type: " + type);
          ASSERT(false);
        }

      }
    }

    loaded_half = working_half;
  }


  // return the datapoint at the requested index
  return data[idx];

}

// takes index in bits
// required byte aligned bit idx
uint32_t datastream::readInt_ds(uint64_t& bit_idx, uint8_t bytes){
  ASSERT(bytes >= 1 && bytes <= 8);
  ASSERT(bit_idx % 8 == 0);

  // make sure everything is loaded nicely
  (void)this->operator[]( (bit_idx / 8) + bytes - 1 );

  // do a byteswap of the bytes
  uint32_t val = 0;
  for(uint8_t i = 0; i < bytes; ++i) val = (val << 8) | data[(bit_idx / 8) + bytes - i - 1];

  bit_idx += 8 * bytes; 
  bit_idx %= wind_sz * 3 * 8;
  return val; 
}

/* Reads no more than 32k bytes from the datastream at a specific index into
 * some provided datastructure
 *
 */
void datastream::readBytes(void* data_vd, uint64_t& bit_idx, uint32_t bytes){
  ASSERT(bytes <= wind_sz);
  uint8_t* data_ptr = reinterpret_cast<uint8_t*>(data_vd);


  // optionally load the datastream so future accesses can be faster 
  // should split off a load function to do only that
  (void)this->operator[]( (bit_idx / 8) + bytes - 1 );
  uint32_t temp_idx;

  // cant use std::copy for circular buffer?
  for(uint32_t i = 0; i < bytes; i++){
    temp_idx = ((bit_idx / 8) + i) % (3 * wind_sz);

    *data_ptr = data[temp_idx];
    data_ptr++;
  }

  bit_idx += 8 * bytes;
  bit_idx %= wind_sz * 3 * 8;
}

/* read 0-8 bits from the datastream at some index
 * aligns the bits to the right side of the char
 * reverses the returned order of the bits in the byte because of deflace LSB
 * see https://stackoverflow.com/questions/2602823/in-c-c-whats-the-simplest-way-to-reverse-the-order-of-bits-in-a-byte for explanation
 * UPDATE: added a boolean to reverse bytes optionally
 */
uint8_t datastream::readBits(uint64_t& bit_idx, uint8_t bits, bool reverse){
  ASSERT(bits <= 8);

  uint32_t byte_idx = (bit_idx / 8);
  uint8_t bit_off = bit_idx % 8;


  uint8_t comp_op = (static_cast<uint16_t>(1) << (8 - bit_off)) - 1;

  // read the byte and reverse bit order
  uint8_t first_byte = this->operator[]( byte_idx );
  first_byte = (first_byte * 0x0202020202ULL & 0x010884422010ULL) % 1023;

  uint8_t second_byte = 0; // could very well leave this as garbage, will be & with 0 if unused
  if(bit_off + bits > 8){
    second_byte = this->operator[]( byte_idx + 1 );
    second_byte = (second_byte * 0x0202020202ULL & 0x010884422010ULL) % 1023;
  }

  uint8_t left_aligned_byte = ((first_byte & comp_op) << bit_off) | ((second_byte & ~comp_op) >> (8 - bit_off));

  bit_idx += bits;
  bit_idx %= wind_sz * 3 * 8;

  left_aligned_byte = left_aligned_byte >> (8 - bits);
  if(reverse) {
    left_aligned_byte = (left_aligned_byte * 0x0202020202ULL & 0x010884422010ULL) % 1023;
    left_aligned_byte = left_aligned_byte >> (8 - bits);
  }

  return left_aligned_byte;
}
