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

constexpr double ugconn = 8314.5;  // J/(kmol-K)
constexpr double light_speed = 2.99792458e8; // m/s
constexpr double pi = 3.14159265358979323846;

constexpr double avogadro = 6.02214076e23; // 1/mol
constexpr double planck = 6.62607015e-34; // m^2 kg / s 
constexpr double boltzmann = 1.380649e-23; // m^2 kg / (s^2 - K)

#define NOW std::chrono::high_resolution_clock::now();

#endif



