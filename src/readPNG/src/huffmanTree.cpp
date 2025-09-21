#include "../include/huffmanTree.hpp"
#include "../include/safeRead.h"
#include <bitset>
#include <memory>
#include <vector>

huffmanTree::huffmanTree():initialized(false){
  codeVals = std::vector<codeType>(MAX_BITS+1);
}

huffmanTree::huffmanTree(int* codeLengths, int length):initialized(true){
  codeVals = std::vector<codeType>(MAX_BITS+1);
  initialize(codeLengths, length);
}

#if 1
void huffmanTree::initialize(int* codeLengths, int length){
  initialized = true;
  min_length = MAX_BITS + 1; // speed up first checks 
  max_length = 0; // speed up failing codes

  int bl_count[MAX_BITS+1] = {0};
  int next_code[MAX_BITS+1] = {0};

  // count the number of codes of each length
  for (int i = 0; i < length; ++i) {
    assert(codeLengths[i] <= MAX_BITS);
    std::cout << codeLengths[i] << ", ";
    bl_count[codeLengths[i]]++;
    if(codeLengths[i] < min_length) min_length = codeLengths[i];
    if(codeLengths[i] > max_length) max_length = codeLengths[i];
  }
  std::cout << std::endl;

  // find the numerical value of the smallest code for each code length
  int code = 0;
  // bl_count[0] = 0; // already 0
  for (int bits = 1; bits <= MAX_BITS; bits++) {
    code = (code + bl_count[bits-1]) << 1;
    codeVals.at(bits).minValue = code;
  }

  for(int i = 0; i < length; i++){
    if(codeLengths[i] == 0) continue;
    assert(codeLengths[i] <= MAX_BITS);
    codeVals.at(codeLengths[i]).values.push_back(i); // improve this later by initializing 
                                                     // vector first and then having a 
                                                     // increment operator to assign
  }
}

uint16_t huffmanTree::decode(datastream& zlib_ds, uint64_t& bit_idx){
  assert(initialized);
  uint16_t bitValue; // PNG files wont have huffman codes > 16 bits
  uint8_t length = std::min(min_length, static_cast<uint8_t>(8)); // cant read more than 8 bits at a time
  bitValue = zlib_ds.readBits(bit_idx, length);

  while(1){
    if(bitValue < codeVals.at(length).minValue); // just to make next if statement more readable
    else if(bitValue - codeVals.at(length).minValue < \
        codeVals.at(length).values.size()) break;

    bitValue = (bitValue << 1) | zlib_ds.readBits(bit_idx, 1);
    length++;
    if(length > max_length) {
      ERROR("invalid code");
      assert(false);
    }
  }

  uint16_t value = codeVals.at(length)[bitValue - codeVals.at(length).minValue];
  // disp(value);
  return value;
}

#else

// not using vector for the codeLengths because I want to support static arrays as well
void huffmanTree::initialize(int* codeLengths, int length){
  initialized = true;
  std::vector<int> huffmanCodes = std::vector<int>(length, -1);

  int bl_count[MAX_BITS+1] = {0};
  int next_code[MAX_BITS+1] = {0};

  // count the number of codes of each length
  for (int i = 0; i < length; ++i) {
    assert(codeLengths[i] <= MAX_BITS);
    bl_count[codeLengths[i]]++;
  }

  // find the numerical value of the smallest code for each code length
  int code = 0;
  // bl_count[0] = 0; // already 0
  for (int bits = 1; bits <= MAX_BITS; bits++) {
    code = (code + bl_count[bits-1]) << 1;
    next_code[bits] = code;
  }

  // generate the huffman codes
  for (int n = 0;  n < length; n++) {
    int len = codeLengths[n];
    if (len != 0) {
      huffmanCodes[n] = next_code[len];
      next_code[len]++;
    }
  }

  // generate a tree based on the codes
  for (int i = 0; i < length; ++i) {
    if (huffmanCodes[i] != -1) {
      this->insertAtBinaryPath(i, huffmanCodes[i], codeLengths[i]);
    }
  }
}

uint16_t huffmanTree::decode(datastream& zlib_ds, uint64_t& bit_idx){
  ASSERT(initialized);
  std::shared_ptr<Node> iter = root;
  uint8_t nextBit;

  while(1){
    // if its a leaf node, we have successfully traversed the tree
    if(!iter->right && !iter->left){
      disp(iter->value);
      return iter->value;
    }

    nextBit = zlib_ds.readBits(bit_idx, 1);

    // continue traversing the tree
    if(nextBit && iter->right){
      iter = iter->right;
    }else if(!nextBit && iter->left) [[likely]] {
      iter = iter->left;
    }else{
      ERROR("Invalid code");
      assert(false);
    }
  }
}
#endif

void huffmanTree::printTree(const std::string& prefix, bool isLeft){
  printTree(root.get(), prefix, isLeft);
}

// Recursive helper to print with indentation
void huffmanTree::printTree(Node* node, const std::string& prefix, bool isLeft){
  if (!node) return;

  // Print current node
  std::cout << prefix;
  std::cout << (isLeft ? "├──" : "└──" );
  std::cout << node->value << std::endl;

  // Recurse on children
  if (node->left) {
    printTree(node->left.get(), prefix + (isLeft ? "│   " : "    "), true);
  }
  if(node->right){
    printTree(node->right.get(), prefix + (isLeft ? "│   " : "    "), false);
  }
}

void huffmanTree::insertAtBinaryPath(int value, int binaryCode, int codeLen) {
  if (!root) {
    root = std::make_shared<Node>(-1); // Initialize root if it doesn't exist
  }

  std::shared_ptr<Node> current = root;

  // We will use binary digits of binaryCode to navigate the tree
  if(!codeLen) return;
  int mask = 1 <<  (codeLen - 1);

  while (mask > 1) {
    if (binaryCode & mask) {  // If current bit is 1, go to the right
      if (!current->right) {
        current->right = std::make_shared<Node>(-1); // Create node if not existing
      }
      current = current->right;
    } else {  // If current bit is 0, go to the left
      if (!current->left) {
        current->left = std::make_shared<Node>(-1); // Create node if not existing
      }
      current = current->left;
    }
    mask >>= 1; // Move to the next bit
  }

  // Last bit: place the value at the final position
  if (binaryCode & 1) {  // If the last bit is 1, go to the right
    if (!current->right) {
      current->right = std::make_shared<Node>(value); // Create final node if needed
    } else {
      current->right->value = value;  // Replace value if node exists
    }
  } else {  // If the last bit is 0, go to the left
    if (!current->left) {
      current->left = std::make_shared<Node>(value); // Create final node if needed
    } else {
      current->left->value = value;  // Replace value if node exists
    }
  }
}

