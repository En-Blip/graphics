#pragma once

#include "vectors.h"
#include "error.h"
#include <algorithm>
#include <iostream>
#include <numbers>

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double e = 2.71828182845904523536028747135266249775724709369995;


int max (int a, int b);
double max (double a, double b);
int min (int a, int b);

template<typename T>
T sq(T a){
    return a*a;
}

// hashing functions
float2 genHash(int2 bounds, float2 seed);

// math functions
float fract(float x);
float dot(float2 A1, float2 A2);
float dot(float3 A1, float3 A2);
float3 normalize(float3 A1);
float2 normalize(float2 A1);

float dist(int2 A1, int2 A2);
float dist(int3 A1, int3 A2);
double dist(double2 A1, double2 A2);

// blending functions
int linearBlend(float input, int2 inputBounds, int2 outputBounds);
float flinearBlend(float input, float2 inputBounds, float2 outputBounds);
float cubicBlend(float input, float2 inputBounds, float2 outputBounds);

// float cubicBlend(float input, float2 inputBounds, float2 outputBounds);


