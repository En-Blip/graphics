#include "../include/readPNG.h"
#include "../include/datastream.hpp"
#include "../include/huffmanTree.hpp"
#include <algorithm>
#include <bitset>
#include <cstdint>
#include <fstream>
#include <vector>

/* IMPROVEMENTS
 * instead of getting all the chunks at once, we can use a 32768 byte sliding window
 * (or whatever is defined) and discard datastream values further than the sliding window back
 * This could be done by making the datastream a circular buffer
 *
 *
 */


// see reference in https://www.rfc-editor.org/rfc/rfc1950,
// https://www.rfc-editor.org/rfc/rfc1951#section-2,
// or https://stackoverflow.com/questions/73779035/how-does-zlib-decompression-work-on-a-png-idat-chunk for the summary
std::vector<color> loadImgFromPng(const std::string& filename, int2& dims){
    const std::vector<color> nullvec; // some vector to return upon failure

    // load image and return the first eight bytes
    // open file
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        ERROR("Error opening the file: " + filename);
        return nullvec;
    }

    // read header
    if(!readPNGheader(file)){
        ERROR("PNG header mismatch");
        return nullvec;
    }

    // find and read IHDR chunk to check settings
    IDHR_t settings;
    if(!readIHDR(file, settings)){
        ERROR("No IHDR chunk");
        return nullvec;
    }

    dims.x = settings.imgDims.x;
    dims.y = settings.imgDims.y;

    // read the CRC (and ignore it for now)
    uint32_t CRC = readInt(file, 4);

    ERROR("currently pretends there is no PLTE chunk");

    // create the datastream object to begin decoding the data
    datastream zlib_ds = datastream(file);
    // create some iterator to keep track of working file idx
    uint64_t bit_idx = 0;
    ERROR("this bit_idx may end up being too small, should have some check to prevent overflow");


    // read zlib header from input stream
    uint8_t CMF = zlib_ds[0];
    uint8_t CM = CMF & 0x0F; // compression method
    uint8_t CINFO = (CMF & 0xF0) >> 4; // compression info 

    ASSERT(CM == 8); // CM should be 8 for PNG files
    ASSERT(CINFO < 8); // values above 7 not allows in this version of the specification

    uint32_t wind_sz = 1u << (CINFO + 8); // or is it pow(2, cinfo) + 8; ??
                                          // window size doesnt really matter for us because I
                                          // want to limit file reads anyways, but may be a 
                                          // good redundancy check in the future

    uint8_t FLG = zlib_ds[1];
    uint8_t FLEVEL = (FLG & 0xc0) >> 6;
    uint8_t FDICT = (FLG & 0x20) >> 5;
    uint8_t FCHECK = FLG & 0x1F;

    assert(((CMF << 8) + FLG) % 31 == 0); // check for validity of CMF and FLG;

    // read dict if required
    uint32_t DICT;
    bit_idx = 2*8;
    if(FDICT){
        DICT = zlib_ds.readInt_ds(bit_idx);
    }

    // decoding specific variables
    bool lastBlock = false;
    std::unique_ptr<uint8_t[]> output_ds = std::make_unique<uint8_t[]>((settings.imgDims.x * 3 + 1) * settings.imgDims.y);
    uint8_t* output_cursor = output_ds.get();


    // read and decode compressed data
    do{
        // read block header from input stream
        uint8_t HDR = zlib_ds.readBits(bit_idx, 3);
        const uint8_t BFINAL = (HDR & 0b100) >> 2;
        const uint8_t BTYPE = HDR & 0b11;

        assert(BTYPE != 0b11);

        lastBlock = BFINAL;

        // if stored with no compression
        if(BTYPE == 0b00){
            // skip any remaining bits in current partially processed byte
            bit_idx += (bit_idx % 8) ? (8 - bit_idx % 8) : 0;

            // read LEN and NLEN
            uint16_t LEN = zlib_ds.readInt_ds(bit_idx, 2);
            uint16_t NLEN = zlib_ds.readInt_ds(bit_idx, 2);

            ASSERT((uint16_t)(~LEN) == NLEN);

            // copy LEN bytes of data to output
            while(LEN > 0){
                uint32_t bytes_read = std::min(LEN, static_cast<uint16_t>(zlib_ds.wind_sz));
                zlib_ds.readBytes(output_cursor, bit_idx, \
                        bytes_read); // LEN will never be more 

                assert((output_cursor - output_ds.get()) < settings.imgDims.x * settings.imgDims.y * 3);

                LEN -= bytes_read;
                output_cursor += bytes_read;
            }
        }else{
            huffmanTree huffmanCodes;
            huffmanTree distTree;

            uint8_t bitval, extraBits;
            uint32_t dist;
            uint32_t length;
            uint32_t value;

            // if compressed with dynamic Huffman codes
            if(BTYPE == 0b10){

                // read the various lengths of the block values
                uint8_t HLIT = zlib_ds.readBits(bit_idx, 5) + 257;
                uint8_t HDIST = zlib_ds.readBits(bit_idx, 5) + 1;
                uint8_t HCLEN = zlib_ds.readBits(bit_idx, 4) + 4;

                uint8_t code_idx[] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
                int codeCodeLens[19] = { 0 };

                // create the huffman tree to decode huffman codes (huffman code inception)
                for(uint8_t i = 0; i < HCLEN; i++) 
                    codeCodeLens[code_idx[i]] = zlib_ds.readBits(bit_idx, 3);

                if(HDIST == 2 && !codeCodeLens[0]) ERROR("Check functionality of special case ");

                huffmanTree codeCodeTree(codeCodeLens, 19);
                codeCodeTree.printTree();

                int previousCodeLength = -1;
                auto codeLens = std::make_unique<int[]>( HLIT );
                auto distLens = std::make_unique<int[]>( HDIST );

                uint8_t repeatNum;
                for (uint16_t i = 0; i < (HLIT);) {
                    value = codeCodeTree.decode(zlib_ds, bit_idx);
                    
                    if(value < 16){
                        previousCodeLength = value;
                        codeLens[i] = previousCodeLength;
                        i++;
                        continue;
                    }else if(value == 16){
                        assert(previousCodeLength >= 0);
                        repeatNum = 3 + zlib_ds.readBits(bit_idx, 2);
                        for(uint8_t x = 0; x < repeatNum; x++){
                            codeLens[i] = previousCodeLength;
                            i++;
                        }
                        continue;
                    }else if(value == 17){
                        repeatNum = 3 + zlib_ds.readBits(bit_idx, 3);
                    }else if(value == 18){ 
                        ERROR("add asserts here to make sure no buffer overflow");
                        repeatNum = 11 + zlib_ds.readBits(bit_idx, 7);
                    }
                    for(uint8_t x = 0; x < repeatNum; x++){
                        codeLens[i] = 0;
                        i++;
                    }
                }

                for (uint16_t i = 0; i < (HDIST);) {
                    value = codeCodeTree.decode(zlib_ds, bit_idx);
                    
                    if(value < 16){
                        previousCodeLength = value;
                        distLens[i] = value;
                        i++;
                    }else if(value == 16){
                        uint8_t repeatNum = 3 + zlib_ds.readBits(bit_idx, 2);
                        for(uint8_t x = 0; x < repeatNum; x++){
                            distLens[i] = previousCodeLength;
                            i++;
                        }
                    }else if(value == 17){
                        uint8_t repeatNum = 3 + zlib_ds.readBits(bit_idx, 3);
                        for(uint8_t x = 0; x < repeatNum; x++){
                            distLens[i] = 0;
                            i++;
                        }
                    }else if(value == 18){
                        uint8_t repeatNum = 11 + zlib_ds.readBits(bit_idx, 7);
                        for(uint8_t x = 0; x < repeatNum; x++){
                            distLens[i] = 0;
                            i++;
                        }
                    }
                }

                huffmanCodes.initialize(codeLens.get(), HLIT + 257);
                distTree.initialize(distLens.get(), HDIST + 1);

            }else{
                int codeLens[] = {8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8};

                // int codeLens[288]; // same thing as above but more readable
                // for(uint16_t i = 0; i < 144; i++) codeLens[i] = 8;
                // for(uint16_t i = 144; i < 256; i++) codeLens[i] = 9;
                // for(uint16_t i = 256; i < 280; i++) codeLens[i] = 7;
                // for(uint16_t i = 280; i < 288; i++) codeLens[i] = 8;

                int distLens[] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

                huffmanCodes.initialize(codeLens, sizeof(codeLens) / sizeof(codeLens[0]));
                distTree.initialize(distLens, sizeof(distLens) / sizeof(distLens[0]));
            }

            while(1){
                value = huffmanCodes.decode(zlib_ds, bit_idx);

                ASSERT(value <= 285); // max value range

                if(value < 256){
                    // copy value (literal byte) to output stream
                    *output_cursor = static_cast<uint8_t>(value);
                    output_cursor++;
                    assert((output_cursor - output_ds.get()) < settings.imgDims.x * settings.imgDims.y * 3);
                }else{
                    if(value == 256) break; // end of block
                                            
                    // use value and extra bits to decode true length
                    if(value == 285) length = 258;
                    else if (value < 265) length = value - 254;
                    else{
                        uint8_t minLen[] = {3, 4, 5, 6, 7, 8, 9 , 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227};
                        extraBits = (value - 261) / 4;

                        length = minLen[value-257] + zlib_ds.readBits(bit_idx, extraBits);
                    }

                    disp(int(length));
                                            
                    // decode distance from input stream
                    uint8_t distCode = distTree.decode(zlib_ds, bit_idx);

                    // use dist value and extra bits to decode the true distance
                    if(distCode >= 4){
                        uint32_t minDist[] = {1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577}; // the minimum distance for each code
                        extraBits = (distCode - 2) / 2;

                        dist = 0;
                        if(extraBits > 8){
                            dist = zlib_ds.readBits(bit_idx, extraBits - 8);
                            dist <<= 8;
                            extraBits -= 8;
                        }
                        dist |= static_cast<uint32_t>(zlib_ds.readBits(bit_idx, extraBits));
                        dist += minDist[distCode];
                    }else{
                        dist = distCode + 1;
                    }

                    // move backwards distance bytes in the output
                    // stream, and copy length bytes from this
                    // position to the output stream.
                    disp(int(distCode));
                    disp(dist);
                    disp(output_cursor - output_ds.get());
                    uint8_t* cpy_cursor = output_cursor - dist;

                    assert(cpy_cursor >= output_ds.get());
                    assert((output_cursor + length - output_ds.get()) < settings.imgDims.x * settings.imgDims.y * 3);

                    memmove(output_cursor, cpy_cursor, length);
                    output_cursor += length;
                }
            }
        }

    }while(!lastBlock);

    // read ADLER32 checksum

    // unfilter data
    // currently we assert no filtering
    
    auto colorBuffer = std::vector<color>(settings.imgDims.x * settings.imgDims.y);

    output_cursor = output_ds.get();
    color tmp_col;
    disp(settings.imgDims.x);
    disp(settings.imgDims.y);
    for(int i = 0; i < settings.imgDims.y; i++){
        disp(int(*output_cursor));
        output_cursor++;

        for(int j = 0; j < settings.imgDims.x; j++){
            tmp_col.r = *output_cursor;
            output_cursor++;
            tmp_col.g = *output_cursor;
            output_cursor++;
            tmp_col.b = *output_cursor;
            output_cursor++;
            colorBuffer.at(i*settings.imgDims.x + j) = tmp_col;
        }
    }

    // create the color buffer to return
    return colorBuffer;
}

bool readPNGheader(std::ifstream& file){
    std::unique_ptr<char[]> header = std::make_unique<char[]>(8);
    char PNGheader[8] = { -119, 80, 78, 71, 13, 10, 26, 10 }; // -119 is 137 as a char
    safeRead(file, header.get(), 8);
    for (int i = 0; i < 8; ++i) {
        if (header[i] != PNGheader[i]) {
            std::cout << header[i] << " != " << PNGheader[i] << std::endl;
            return false;
        }
    }
    return true;
}

bool readIHDR(std::ifstream& file, IDHR_t& settings){
    uint32_t chunkSize;
    safeRead(file, &chunkSize, 4);

    // read the chunk type, since its the first chunk it should be IHDR
    char chunkType[4];
    safeRead(file, chunkType, 4);
    if (chunkType[0] != 'I' || chunkType[1] != 'H' || chunkType[2] != 'D' || chunkType[3] != 'R') {
        return false;
    }

    settings.imgDims.x = readInt(file, 4);
    settings.imgDims.y = readInt(file, 4); 

    // if this is IHDR, will contain, in bytes
    // width(4), height(4), bit depth(1), color type(1), compression method(1)
    // filter method(1), interlace method(1)
    safeRead(file, &settings.bit_depth, 5); // careful: cpp compiler inserts more bited in struct
                                            // for alignment -> dont use sizeof

                                            // make sure the settings work for what I have 
    // ASSERT(settings.imgDims.x * settings.imgDims.y == 3024 * 1964);
    // ASSERT(settings.color_type == 2);
    ASSERT(settings.bit_depth == 8 || settings.bit_depth == 16);
    ASSERT(settings.compression_method == 0);
    ASSERT(settings.filter_method == 0);
    ASSERT(settings.interlace_method == 0);

    return true;
}
