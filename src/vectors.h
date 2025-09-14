#pragma once
#include <cstdint>

struct int2{
    int x;
    int y;
};
struct alignas(4) uint2{uint32_t x; uint32_t y;};
struct alignas(4) float2{float x; float y;};
struct alignas(8) double2{double x; double y;};
struct alignas(4) float3{float x; float y; float z;};
struct alignas(4) int3{int x; int y; int z;};
struct alignas(1) color{unsigned char r; unsigned char g; unsigned char b;};
