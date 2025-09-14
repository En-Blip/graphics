#include "image.h"
#include "readPNG/include/safeRead.h"
#include <bitset>
#include <ios>
#include <iostream>
#include <memory>
#include <vector>

void Tree::insertAtBinaryPath(int value, int binaryCode, int codeLen) {
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

// Recursive helper to print with indentation
void Tree::printTree(std::shared_ptr<Node> node, const std::string& prefix, bool isLeft) {
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
 
void Tree::print(){
    std::bitset<MAX_BITS+1> bit_seq;
    if(root) print(root, bit_seq, 0);
}

void Tree::print(std::shared_ptr<Node> node, std::bitset<MAX_BITS+1>& bit_seq, uint8_t cur_depth){
    assert(cur_depth <= MAX_BITS + 1);

    if(!node->left && !node->right){
        for(int i = cur_depth; i < MAX_BITS+1; i++) bit_seq.reset(i);
        bit_seq.set(cur_depth);
        disp(bit_seq);
        disp(node->value);
        std::cout << "\n" << std::endl;
        return;
    }

    if(node->left){
        bit_seq.reset(cur_depth);
        print(node->left, bit_seq, cur_depth+1);
    }    
    if(node->right){
        bit_seq.set(cur_depth);
        print(node->right, bit_seq, cur_depth+1);
    }


}

int getBitRange(unsigned char* bytes, int bitStart, int bitEnd, int byteSize) {
    if(bitStart > bitEnd || bitEnd > byteSize * 8 || bitStart < 0) { 
        std::cerr << "Bit range error" << std::endl;
        return 0;
    }

    int bits = 0;
    for (int i = bitStart; i <= bitEnd; ++i) {
        bits = (bits << 1) | ((bytes[i / 8] >> (i % 8)) & 1);
    }
    return bits;
}

// reads up to 32 bits and returns their value, then increments the cursor
unsigned int readBits(unsigned char* bytes, int& bitIdx, int numBits, int byteSize){
    unsigned char bits = (unsigned char)getBitRange(bytes, bitIdx, bitIdx + numBits, byteSize);
    bitIdx += numBits;
    return bits;
}

Tree* genHuffmanCodes(int* codeLengths, int length){ // requires code lengths to be alphabet size
    ERROR("refactor this function, this causes a memory leak in the way we're using it atm");
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
    Tree* tree = new Tree();


    // generate a tree based on the codes
    for (int i = 0; i < length; ++i) {
        if (huffmanCodes[i] != -1) {
            // std::cout << i << " " << std::bitset<8>(huffmanCodes[i]) << " " << codeLengths[i] << std::endl;
            tree->insertAtBinaryPath(i, huffmanCodes[i], codeLengths[i]);
        }
    }
    return tree;
}

void writeImgToCSV(const int* data, const std::string& filename, const int2 dimensions) {
    std::ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return;
    }

    for (int i = 0; i < dimensions.y; ++i) {
        for (int j = 0; j < dimensions.x; ++j) {
            outputFile << data[dimensions.x * i + j];

            if (j != dimensions.x - 1) {
                outputFile << ",";
            }
        }
        outputFile << "\n";
    }

    outputFile.close();
}

/**
 * @brief Loads an image from a PNG file
 * @note PNG reference: http://www.libpng.org/pub/png/spec/1.2/PNG-Chunks.html
 * @note only supported color type as of now is RGB
 * @param filename The filename to load
 * @return The image data in a color array
 */
color* loadImgFromPng1(const std::string& filename) {
    // load image and return the first eight bytes

    // open file
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return nullptr;
    }

    // read header
    std::unique_ptr<char[]> header = std::make_unique<char[]>(8);
    char PNGheader[8] = { -119, 80, 78, 71, 13, 10, 26, 10 }; // -119 is 137 as a char
    file.read(header.get(), 8);
    for (int i = 0; i < 8; ++i) {
        if (header[i] != PNGheader[i]) {
            std::cout << header[i] << " != " << PNGheader[i] << std::endl;
            ERROR("PNG header mismatch");
            return nullptr;
        }
    }

    // find IHDR chunk
    // read the length bytes
    uint32_t chunkSize = readInt(file, 4);

    // read the chunk type, since its the first chunk it should be IHDR
    char chunkType[4];
    file.read(chunkType, 4);
    if (chunkType[0] != 'I' || chunkType[1] != 'H' || chunkType[2] != 'D' || chunkType[3] != 'R') {
        ERROR("No IHDR chunk");
        return nullptr;
    }

    // if this is IHDR, will contain, in bytes
    // width(4), height(4), bit depth(1), color type(1), compression method(1)
    // filter method(1), interlace method(1)
    int2 imgDims = {readInt(file, 4), readInt(file, 4)};
    char bit_depth, color_type, compression_method, filter_method, interlace_method;
    file.read(&bit_depth, 1); // theres probably a better way to do this
    file.read(&color_type, 1);
    file.read(&compression_method, 1);
    file.read(&filter_method, 1);
    file.read(&interlace_method, 1);

    // make sure the settings work for what I have 
    ASSERT(imgDims.x * imgDims.y == 3024 * 1964);
    ASSERT(color_type == 2);
    ASSERT(bit_depth == 8 || bit_depth == 16);
    ASSERT(compression_method == 0);
    ASSERT(filter_method == 0);
    ASSERT(interlace_method == 0);

    uint32_t CRC = readInt(file, 4);

    // create the color buffer
    color* colorBuffer = new color[imgDims.x * imgDims.y];
    color* cBuffCpy = colorBuffer; // for copying data into the buffer

    unsigned char* datastream = new unsigned char[imgDims.x * imgDims.y * 3];
    unsigned char* datastream_ptr = datastream;

    unsigned char* outputstream = new unsigned char[imgDims.x * imgDims.y * 3];
    unsigned char* outputstream_ptr = outputstream;

    while(1){
        // loop through chunks until IEND
        chunkSize = readInt(file, 4);
        char* chunkType = new char[4];
        file.read(chunkType, 4);

        // check the chunk type
        if(chunkType[0] == 'I' && chunkType[1] == 'E' && chunkType[2] == 'N' && chunkType[3] == 'D'){
            break;
        }else if(chunkType[0] == 'I' && chunkType[1] == 'D' && chunkType[2] == 'A' && chunkType[3] == 'T'){
            // IDAT chunk is correct so do nothing
        }else{
            std::cerr << "Unknown chunk type: " << chunkType << std::endl;
            return colorBuffer;
        }

        delete[] chunkType;

        // read the data
        char* dataBuffer = new char[chunkSize];
        file.read(dataBuffer, chunkSize);

        // concatenate the data to the zlib data-stream
        for (int i = 0; i < chunkSize; ++i) {
            *datastream_ptr = dataBuffer[i];
            ++datastream_ptr;
        }

        delete[] dataBuffer;
        CRC = readInt(file, 4);
    }

    // when reached IEND chunk return
    std::cout << "Reached IEND chunk, decoding ..." << std::endl;

    // close file
    file.close();

    // decompress the data
    int sizeofDatastream = datastream_ptr - datastream;

    // check if the compression method is 8 (deflate/inflate)
    char zlib_comp_method = datastream[0] & 0x0F; // first 4 bits are method
    char compInfo = datastream[0] & 0xF0 >> 4; // next 4 bits are further info
    int windowSize = pow(2, compInfo) + 8;
    ASSERT(zlib_comp_method == 8);

    /**@note next bit is flag
     * 0-4 FCHECK
     * 5 FDICT
     * 6-7 FLEVEL (see https://www.rfc-editor.org/rfc/rfc1950#section-4)
    */
    char FDICT = datastream[1] & 0b00010000;                                // is this where datasream should be starting (and not from the byte boundary?)

    int bitIdx = 0;
    unsigned char* data = datastream + 2;
    char BFINAL;
    char BTYPE;
    do{
        // first 3 bits are the header
        BFINAL = (char)readBits(data, bitIdx, 1, sizeofDatastream);
        BTYPE = (char)readBits(data, bitIdx, 2, sizeofDatastream);

        if(!BTYPE) { // stored without compression
            // skip any unprocessed bits
            readBits(data, bitIdx, 5, sizeofDatastream);
            // read LEN and NLEN
            int LEN = readBits(data, bitIdx, 8, sizeofDatastream) + 256 * readBits(data, bitIdx, 8, sizeofDatastream); // this was originally switched *********** (like the second byte was shifted -> little endian)
            int NLEN = readBits(data, bitIdx, 8, sizeofDatastream) + 256 * readBits(data, bitIdx, 8, sizeofDatastream);
            ASSERT(LEN == (NLEN ^ 0xFFFF));

            // copy LEN bytes
            for(int j = 0; j < LEN; j++){
                *outputstream = readBits(data, bitIdx, 8, sizeofDatastream);
                outputstream++;
            }
        }else{
            Tree* huffmanCodes = nullptr;
            if(BTYPE == 2) { // dynamic Huffman codes
                // read representation of code trees
                // SEE SUBSECTION BELOW
                ERROR("image.cpp: dynamic duffman codes are not implemented"); // ^^
                return nullptr; 
            }else{
                int codeLens[] = {8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8};
                huffmanCodes = genHuffmanCodes(codeLens, ALPHABET_SIZE);
            }
            unsigned int value = 0;
            char bitval = 0, extraBits = 0;
            unsigned int dist = 0;
            unsigned int length = 0;
            std::shared_ptr<Node> root = huffmanCodes->root;
            std::shared_ptr<Node> current = root;

            int distLens[] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
            std::shared_ptr<Node> distRoot = genHuffmanCodes(distLens, 29)->root;
            while(value != 256){
                // decode literal length/value from dstream
                ASSERT(current);
                while(1){
                    bitval = readBits(data, bitIdx, 1, sizeofDatastream);
                    if(bitval == 0){
                        current = current->left;
                    }else{
                        current = current->right;
                    }
                    if(!current->left && !current->right){
                        value = current->value;
                        break;
                    }
                }

                // if value < 256, copy to output stream
                if(value < 256){
                    *outputstream = value;
                    outputstream++;
                }else{
                    if(value == 256){
                        break;
                    }else{
                        //decode length from input stream
                        if (value >= 257 && value <= 264) {
                            length = 3 + (value - 257);
                            extraBits = 0;
                        } else if (value >= 265 && value <= 268) {
                            length = 11 + 2 * (value - 265);
                            extraBits = 1;
                        } else if (value >= 269 && value <= 272) {
                            length = 19 + 4 * (value - 269);
                            extraBits = 2;
                        } else if (value >= 273 && value <= 276) {
                            length = 35 + 8 * (value - 273);
                            extraBits = 3;
                        } else if (value >= 277 && value <= 280) {
                            length = 67 + 16 * (value - 277);
                            extraBits = 4;
                        } else if (value >= 281 && value <= 284) {
                            length = 131 + 32 * (value - 281);
                            extraBits = 5;
                        } else if (value == 285) {
                            length = 258;
                            extraBits = 0;
                        } else {
                            // Invalid length code
                            ERROR("Invalid length code");
                        }
                        length += readBits(data, bitIdx, extraBits, sizeofDatastream);

                        current = distRoot;
                        ASSERT(current);
                        while(1){
                            bitval = readBits(data, bitIdx, 1, sizeofDatastream);
                            if(bitval == 0){
                                current = current->left;
                            }else{
                                current = current->right;
                            }
                            if(!current->left && !current->right){
                                dist = current->value;
                                break;
                            }
                        }

                        if (dist >= 0 && dist <= 3) {
                            dist = 1 + dist;
                            extraBits = 0;
                        } else if (dist >= 4 && dist <= 5) {
                            dist = 5 + 2 * (dist - 4);
                            extraBits = 1;
                        } else if (dist >= 6 && dist <= 7) {
                            dist = 9 + 4 * (dist - 6);
                            extraBits = 2;
                        } else if (dist >= 8 && dist <= 9) {
                            dist = 17 + 8 * (dist - 8);
                            extraBits = 3;
                        } else if (dist >= 10 && dist <= 11) {
                            dist = 33 + 16 * (dist - 10);
                            extraBits = 4;
                        } else if (dist >= 12 && dist <= 13) {
                            dist = 65 + 32 * (dist - 12);
                            extraBits = 5;
                        } else if (dist >= 14 && dist <= 15) {
                            dist = 129 + 64 * (dist - 14);
                            extraBits = 6;
                        } else if (dist >= 16 && dist <= 17) {
                            dist = 257 + 128 * (dist - 16);
                            extraBits = 7;
                        } else if (dist >= 18 && dist <= 19) {
                            dist = 513 + 256 * (dist - 18);
                            extraBits = 8;
                        } else if (dist >= 20 && dist <= 21) {
                            dist = 1025 + 512 * (dist - 20);
                            extraBits = 9;
                        } else if (dist >= 22 && dist <= 23) {
                            dist = 2049 + 1024 * (dist - 22);
                            extraBits = 10;
                        } else if (dist >= 24 && dist <= 25) {
                            dist = 4097 + 2048 * (dist - 24);
                            extraBits = 11;
                        } else if (dist >= 26 && dist <= 27) {
                            dist = 8193 + 4096 * (dist - 26);
                            extraBits = 12;
                        } else if (dist >= 28 && dist <= 29) {
                            dist = 16385 + 8192 * (dist - 28);
                            extraBits = 13;
                        } else {
                            // Invalid distance code
                            ERROR("Invalid distance code");
                        }

                        // Read the extra bits and add them to the base distance
                        int extra = (extraBits > 0) ? readBits(data, bitIdx, extraBits, sizeofDatastream) : 0;
                        dist += extra;
                    }

                    for (int i = dist; i > 0; i--) {
                        *outputstream = *(outputstream - i);
                        outputstream++;
                    }
                    current = root; // is this it
                }
            }
        }
    }while(BFINAL == 0);
    /*
    do
        read block header from input stream.
        if stored with no compression
            skip any remaining bits in current partially
                processed byte
            read LEN and NLEN (see next section)
            copy LEN bytes of data to output
        otherwise
            if compressed with dynamic Huffman codes
                read representation of code trees (see
                subsection below)
            loop (until end of block code recognized)
                decode literal/length value from input stream
                if value < 256
                copy value (literal byte) to output stream
                otherwise
                if value = end of block (256)
                    break from loop
                otherwise (value = 257..285)
                    decode distance from input stream

                    move backwards distance bytes in the output
                    stream, and copy length bytes from this
                    position to the output stream.
            end loop
    while not last block
    */
    return colorBuffer;
}

void convertPNGToChar(const std::string& filein, const std::string& fileout) {
    // open file
    std::ifstream file(filein, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening the file: " << filein << std::endl;
        return;
    }

    std::ofstream outputFile(fileout, std::ios::binary);
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file: " << fileout << std::endl;
        return;
    }
    char to_write;

    int i = 0;
    while(!file.eof() && i < 100){
        // write byte of input file to output file as a character
        to_write = (char)file.get();
        if(to_write == 'I' || to_write == 'H' || to_write == 'D' || to_write == 'R' || to_write == 'E' || to_write == 'N' || to_write == 'A' || to_write == 'T'){
            outputFile << to_write << ",";
        }else{
            outputFile << (int)to_write << ",";
        }

        i++;
    }

    file.close();
    outputFile.close();
}

