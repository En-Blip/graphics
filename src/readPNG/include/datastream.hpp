#pragma once
#include <fstream>
#include "../../error.h"

// #define CHUNKS
// #define VERBOSE

class datastream {
  public: 
    static constexpr uint32_t wind_sz = 32768;

    datastream(std::ifstream& filepath):file(filepath), loaded_half(-1), \
                                        bytes_left_in_chunk(0), end_idx(-1){
                                          assert(file.is_open());
                                        };
    uint8_t operator[](int idx); // used as a loading function as well
    uint32_t readInt_ds(uint64_t& bit_idx, uint8_t bytes = 4);
    void readBytes(void* data_vd, uint64_t& bit_idx, uint32_t bytes);
    void readBits(void* data_vd, uint32_t bit_idx, uint32_t bits);
    uint8_t readBits(uint64_t&  bit_idx, uint8_t bits, bool reverse = false);

  private:
    datastream(); // a file must be provided for this class to work
    datastream(const datastream& src); // prevent copy construction
    datastream& operator=(const datastream& src); // prevent assignment
    datastream(const datastream&& src); // prevent move construction

    std::ifstream& file;
    uint32_t bytes_left_in_chunk; // the datastream reads in sets of 32k bytes, which may not
                                  // align with a chunk boundary. Keep track of this
    int8_t loaded_half; // the half of the datastream filled with the working data
                        // the other half is jsut back reference
                        // -1 is no data, 0 is first third (3rd is max window), 1 is second
                        // third (1st is max window), 2 is third third (2nd is max window)

    int32_t end_idx;  // if the program reaches IEND, it should throw an error for data
                      // accesses beyond this index

    unsigned char data[wind_sz * 3]; // the data array to manage
};
