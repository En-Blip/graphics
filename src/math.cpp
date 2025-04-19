#include "math.h"

float2 genHash(int2 bounds, float2 seed){
    float2 M1 = {3.1251, 17.8737};
    float M2 = 43758.545312;

    float2 v;
    float2 seed2 = {seed.x+3, seed.y+3};
    v.x = (sin((seed.x*M1.x + seed.y*M1.y) * M2)) - static_cast<int>(sin((seed.x*M1.x + seed.y*M1.y) * M2));
    v.y = (sin((seed2.x*M1.x + seed2.y*M1.y) * M2)) - static_cast<int>(sin((seed2.x*M1.x + seed2.y*M1.y) * M2));

    v.x = v.x * bounds.x/2 + bounds.x/2;
    v.y = v.y * bounds.y/2 + bounds.y/2;

    return v;
}

float fract(float x){
    return x - static_cast<int>(x);
}

float dot(float2 A1, float2 A2){
    return (A1.x * A2.x + A1.y * A2.y);
}

float dot (float3 A1, float3 A2){
    return (A1.x * A2.x + A1.y * A2.y + A1.z * A2.z);
}

float3 normalize(float3 A1){
    float val = sqrt(A1.x * A1.x + A1.y * A1.y + A1.z * A1.z);
    return {A1.x/val, A1.y/val, A1.z/val};
}

float2 normalize(float2 A1){  
    float val = sqrt(A1.x * A1.x + A1.y * A1.y);
    return {A1.x/val, A1.y/val};
}

int linearBlend(float input, int2 inputBounds, int2 outputBounds){
    if(inputBounds.x - inputBounds.y == 0){
        return outputBounds.x;
    }
    // create a linear function to map the 2 inputs/outputs
    return static_cast<int>((outputBounds.y-outputBounds.x) * (input - inputBounds.x) / (inputBounds.y - inputBounds.x) + outputBounds.x);
}

float flinearBlend(float input, float2 inputBounds, float2 outputBounds){
    if (input < inputBounds.x || input > inputBounds.y){
        std::cerr << "math.cpp: flinearBlend: input " << input << " out of bounds" << std::endl;
        exit(1);
    }
    if(inputBounds.x - inputBounds.y == 0){
        return outputBounds.x;
    }
    // create a linear function to map the 2 inputs/outputs
    return (outputBounds.y-outputBounds.x) * (input - inputBounds.x) / (inputBounds.y - inputBounds.x) + outputBounds.x;
}

float cubicBlend(float input, float2 inputBounds, float2 outputBounds){
    if(inputBounds.x - inputBounds.y == 0){
        return outputBounds.x;
    }
    // create a cubic function to map the 2 inputs/outputs
    return (outputBounds.y-outputBounds.x) * (-2*pow((input - inputBounds.x) / (inputBounds.y - inputBounds.x), 3) + 3 * pow((input - inputBounds.x) / (inputBounds.y - inputBounds.x), 2)) + outputBounds.x;
}

float dist(int2 A1, int2 A2){
    return sqrt((A1.x - A2.x)*(A1.x - A2.x) + (A1.y - A2.y)*(A1.y - A2.y));
}

float dist(int3 A1, int3 A2){
    return sqrt((A1.x - A2.x)*(A1.x - A2.x) + (A1.y - A2.y)*(A1.y - A2.y) + (A1.z - A2.z)*(A1.z - A2.z));
}

int max (int a, int b){
    return a > b ? a : b;
}

double max (double a, double b){
    return a > b ? a : b;
}