#pragma once

#include "vectors.h"
#include "error.h"
#include <sstream>
#include <fstream>
#include <iostream>

static constexpr uint16_t ALPHABET_SIZE = 287;
static constexpr uint8_t MAX_BITS = 19; //9;

class Node{
  public:
    int value;
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;

    Node(int val):value(val), left(nullptr), right(nullptr){}
};

class Tree{
  public:
    Tree():root(nullptr){}
    void insertAtBinaryPath(int value, int binaryCode, int codeLen);
    void print();
    void print(std::shared_ptr<Node> node, std::bitset<MAX_BITS+1>& bit_seq, \
        uint8_t cur_depth);
    std::shared_ptr<Node> root;
    void printTree(std::shared_ptr<Node> node, const std::string& prefix = "", bool isLeft = true);

  private:
    void print(std::shared_ptr<Node>(node));
    int findBinaryLength(int binaryCode);
};

int getBitRange(unsigned char* bytes, int bitStart, int bitEnd, int byteSize);
unsigned int readBits(unsigned char* bytes, int& bitIdx, int numBits, int byteSize);
Tree* genHuffmanCodes(int* codeLengths, const int length); 
void writeImgToCSV(const int* data, const std::string& filename, const int2 dimensions);
color* loadImgFromPng1(const std::string& filename);
void convertPNGToChar(const std::string& filein, const std::string& fileout);
