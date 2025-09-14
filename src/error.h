#pragma once

#include <iostream>
#include <signal.h>
#include <cassert>

#define ERROR(msg) std::cerr << "Error: " << msg << std::endl
#define ASSERT(x) assert(x)
// #define ASSERT(x) if (!(x)) { \
//     ERROR("Assertion failed: " << #x); \
// }

#define disp(x) std::cout << (#x) << ": " << (x) << std::endl;
