#pragma once

#include <iostream>
#include <signal.h>

#define ERROR(msg) std::cerr << "Error: " << msg << std::endl
#define ASSERT(x) if (!(x)) ERROR("Assertion failed: " << #x)
