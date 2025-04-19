#pragma once

#include "vectors.h"
#include "error.h"
#include <sstream>
#include <fstream>
#include <iostream>

#define ALPHABET_SIZE 287
#define MAX_BITS 9

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
        void insertAtBinaryPath(int value, int binaryCode);
        std::shared_ptr<Node> root;

    private:
        int findBinaryLength(int binaryCode);
};

int readInt(std::ifstream& file, unsigned short int bytes);
int getBitRange(unsigned char* bytes, int bitStart, int bitEnd, int byteSize);
unsigned int readBits(unsigned char* bytes, int& bitIdx, int numBits, int byteSize);
Tree* genHuffmanCodes(int* codeLengths);
void writeImgToCSV(const int* data, const std::string& filename, const int2 dimensions);
color* loadImgFromPng(const std::string& filename);
void convertPNGToChar(const std::string& filein, const std::string& fileout);
