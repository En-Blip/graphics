#include "../include/huffmanTree.hpp"
#include "../include/safeRead.h"
#include <memory>
#include <vector>

huffmanTree::huffmanTree():initialized(false){

}

huffmanTree::huffmanTree(int* codeLengths, int length):initialized(true){
    initialize(codeLengths, length);
}

void huffmanTree::initialize(int* codeLengths, int length){
    initialized = true;
    std::vector<int> huffmanCodes = std::vector<int>(length, 0);

    int bl_count[MAX_BITS+1] = {0};
    int next_code[MAX_BITS+1] = {0};

    // count the number of codes of each length
    for (int i = 0; i < length; ++i) {
        ASSERT(codeLengths[i] <= MAX_BITS);
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
        }else{
            huffmanCodes[n] = -1;
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

void huffmanTree::printTree(const std::string& prefix, bool isLeft){
    printTree(root, prefix, isLeft);
}

// Recursive helper to print with indentation
void huffmanTree::printTree(std::shared_ptr<Node> node, const std::string& prefix, bool isLeft){
    if (!node) return;

    // Print current node
    std::cout << prefix;
    std::cout << (isLeft ? "├──" : "└──" );
    std::cout << node->value << std::endl;

    // Recurse on children
    if (node->left) {
        printTree(node->left, prefix + (isLeft ? "│   " : "    "), true);
    }
    if(node->right){
        printTree(node->right, prefix + (isLeft ? "│   " : "    "), false);
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

