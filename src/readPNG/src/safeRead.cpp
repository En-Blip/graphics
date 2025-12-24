#include "../include/safeRead.h"

bool safeRead(std::ifstream& file, void* data, uint32_t size){
  assert(file.is_open());
  file.read(reinterpret_cast<char*>(data), size); // security issue?
  if (!file) {
    std::cerr << "Failed to read " << size << " bytes from file" << std::endl;
    assert(false);
    // return false;
  }
  return true;
}

int readInt(std::ifstream& file, unsigned short int bytes) {
  ASSERT(file.is_open());

  int val = 0;
  unsigned char byte;
  for (int i = 0; i < bytes; ++i) {
    safeRead(file, &byte, 1);
    val = (val << 8) | byte;
  }
  return val;
}
