#pragma once
#include "datastream.hpp"

constexpr uint8_t MAX_BITS = 9;

class Node{
    public:
        int value;
        std::shared_ptr<Node> left;
        std::shared_ptr<Node> right;

        Node(int val):value(val), left(nullptr), right(nullptr){}
};

class huffmanTree{
    public:
        huffmanTree();
        huffmanTree(int* codeLengths, int length);

        void initialize(int* codeLengths, int length);
        uint16_t decode(datastream& zlib_ds, uint64_t& bit_idx);
        void printTree(std::shared_ptr<Node> node, const std::string& prefix = "", bool isLeft = true);
        void printTree(const std::string& prefix = "", bool isLeft = true);

    private:
        bool initialized; // cant read anything if tree uninitialized
                          // but we may want to declare it in an outer scope
                          // before initialization

        huffmanTree(const huffmanTree& src); // prevent copy construction
        huffmanTree& operator=(const huffmanTree& rhs); // prevent assignment

        void insertAtBinaryPath(int value, int binaryCode, int codeLen);

        // tree variables
        std::shared_ptr<Node> root;
};
