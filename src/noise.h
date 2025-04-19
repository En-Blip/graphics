#pragma once

#include "vectors.h"
#include "math.h"
#include <algorithm>
#include <iostream>

class Noise
{
private:
    int2 dims;
    int layers;
    void (*freqFunction);
    int2 offset;
    float2 vectors[8] = {{1, 1}, {-1,1}, {1,-1}, {-1,-1}, {static_cast<float>(sqrt(2)), 0}, {static_cast<float>(-sqrt(2)), 0}, {0, static_cast<float>(sqrt(2))}, {0, static_cast<float>(-sqrt(2))}};
    float* customConvolveFilter = new float[9];

public:
    Noise(int2 dims);
    ~Noise();

    unsigned char* perlinNoise(int2 dims, int2 cellcount);
    unsigned char* worleyNoise(int2 dims, int frequency);

    void convolve(int2 dims, unsigned char* image, char filter);

    unsigned char* genPerlin(int layers, int2 frequency, double (*weightCalc)(int nLayers, int layer), bool reversed);
    unsigned char* genWorley(int layers, int tileSize, int2 frequency, double (*weightCalc)(int nLayers, int layer), bool reversed);

    unsigned char* gridTransform(unsigned char* domain, int multiplier);
};
