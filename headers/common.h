#ifndef COMMON_H
#define COMMON_H


#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <chrono>

constexpr double bk = 1.380649e-23; // m^2 kg / (s^2 - K)
constexpr double ugconn = 8314.5;  // J/(mol-K)

#define NOW std::chrono::high_resolution_clock::now();

#endif



