#include "../include/datastream.hpp"
#include "../include/safeRead.h"
#include <algorithm>
#include <bitset>
#include <cstddef>
#include "../../image.h"

// TODO: split off a load function from operator[] then use that in readBytes

// DOESNT REVERSE BIT ORDER
uint8_t datastream::operator[](int idx){
  // convert the index into a multiple of the size of the array
  idx %= wind_sz * 3;

  // find which half the access is requesting    
  int8_t working_half = idx / wind_sz;

  // return an error if we are accessing beyond end_idx
  uint8_t end_working_half = end_idx / wind_sz;
  if(end_idx >= 0 && working_half == end_working_half && idx > end_idx){
    ERROR("Access out of datastream bounds");
    ASSERT(working_half != end_working_half && idx <= end_idx);
  }

  // load the appropriate half if needed
  if(loaded_half < working_half || loaded_half == 3 && working_half == 1){ // im gonna ignore
                                                                           // this magic number

                                                                           // keep track of bytes read in datastream
    uint16_t bytes_read = 0;

    while(bytes_read < wind_sz){
      // read the remaining bytes in chunk
      unsigned char* load_point = data + wind_sz * working_half;
      uint32_t bytes_to_read = std::min(bytes_left_in_chunk, wind_sz);
      if(bytes_to_read){
        bytes_read += bytes_to_read;
        bytes_left_in_chunk -= bytes_to_read;
        safeRead(file, load_point, bytes_to_read);
      }

      if(bytes_read < wind_sz){ 
        uint32_t chunkSize = readInt(file, 4);
        disp(chunkSize);

        // read the chunk type
        char chunkType[5] = {0};
        file.read(chunkType, 4);

        disp(chunkType);

        // confirm that it is an IDAT chunk
        if(chunkType[0] == 'I' && chunkType[1] == 'E' && chunkType[2] == 'N' && chunkType[3] == 'D'){
          end_idx = wind_sz * working_half + bytes_read;
          if(working_half == end_working_half && idx > end_idx){
            ERROR("Access out of datastream bounds");
            ASSERT(working_half != end_working_half && idx <= end_idx);
          }
          break;
        }else if(!strcmp(chunkType, "IDAT")){
          // read the data and return it
          bytes_to_read = std::min(chunkSize, wind_sz - bytes_read);
          bytes_left_in_chunk = chunkSize - (bytes_to_read);
          safeRead(file, load_point + bytes_read, bytes_to_read);
          bytes_read += bytes_to_read;

          // read the chunk CRC (and ignore it)
          uint32_t CRC = readInt(file, 4);

        }else if(chunkType[0] & 0x20){
          std::string type(chunkType);
          ERROR("Skipping chunk type: " + type);
          file.seekg(chunkSize + 4, std::ios::cur);
          return this->operator[](idx);
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
  ASSERT(bit_idx / 8);
  // make sure everything is loaded nicely
  (void)this->operator[]( (bit_idx / 8) + bytes - 1 );

  // do a byteswap of the bytes
  uint32_t val = 0;
  for(uint8_t i = 0; i < bytes; ++i) val = (val << 8) | data[(bit_idx / 8) + bytes - i - 1];

  bit_idx += 8 * bytes; 
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
}

/* read n bits into some datastructure aligned with the last byte of that datastructure
 *
 */
void datastream::readBits(void* data_vd, uint32_t bit_idx, uint32_t bits){
  ERROR("function incomplete");
  assert(false);
  ASSERT(bits / 8 < wind_sz); // still dont want to read more than 32 =k bytes

  uint8_t* data_ptr = reinterpret_cast<uint8_t*>(data_vd);
  uint8_t read_bit_off = bits % 8; // compiler should make this equivalent to bit_idx & 0x7;
  uint8_t write_byte = (bit_idx / 8) % (3 * wind_sz);

  // load the database at max dist to speed up next reads
  (void)this[(bit_idx + bits - 1) / 8];


  // read the first not entirely full byte into data
  uint8_t comp_op = (1 << (8 - read_bit_off)) - 1;
  *data_ptr = *data_ptr & (comp_op ^ 1);
  *data_ptr = *data_ptr | (data[write_byte] & comp_op);
  data_ptr++;
  write_byte++;
  bits -= (8 - read_bit_off);


  // read bits 8 at a time and put into data_vd
  while(bits >= 8){
    *data_ptr = data[write_byte];

    data_ptr++;
    write_byte++;
    bits -= 8;
  }

  ASSERT(bits == 0);
}

/* read 0-8 bits from the datastream at some index
 * aligns the bits to the right side of the char
 * reverses the returned order of the bits in the byte because of deflace LSB
 * see https://stackoverflow.com/questions/2602823/in-c-c-whats-the-simplest-way-to-reverse-the-order-of-bits-in-a-byte for explanation
 * UPDATE: added a boolean to reverse bytes optionally
 */
/* validating the algorithm because something works and I cant figure it out for the life of me
 *  | <765>43210 | 765432<10> |
 *  read 5 bits offset of 5
 *
 *  <56701>
 *  
 *  algorithm:
 *
 *  comp_op = 0b1000 - 1;
 *  comp_op = 0b0111;
 *
 *  first_byte = 76543210
 *  second_byte = 76543210
 *
 *  (first_byte & comp_op) << bit_off = 56700000
 *  (second_byte & ~comp_op) << 8-bit_off = 01234
 *
 *  56701
 *  */
uint8_t datastream::readBits(uint64_t& bit_idx, uint8_t bits, bool reverse){
  ASSERT(bits <= 8);

  uint32_t byte_idx = (bit_idx / 8) % (3 * wind_sz);
  uint8_t bit_off = bit_idx % 8;

  
  uint8_t comp_op = (static_cast<uint16_t>(1) << (8 - bit_off)) - 1;

  uint8_t first_byte = this->operator[]( byte_idx );
  uint8_t second_byte = 0; // could very well leave this as garbage, will be & with 0 if unused
  if(bit_off + bits > 8){
    second_byte = this->operator[]( byte_idx + 1 );
    if(reverse) second_byte = (second_byte * 0x0202020202ULL & 0x010884422010ULL) % 1023;
  }

  // reverse bit order
  if(reverse) first_byte = (first_byte * 0x0202020202ULL & 0x010884422010ULL) % 1023;

  uint8_t left_aligned_byte = ((first_byte & comp_op) << bit_off) | ((second_byte & ~comp_op) >> (8 - bit_off));

  bit_idx += bits;
  return left_aligned_byte >> (8 - bits);
}
