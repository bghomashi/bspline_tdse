#pragma once

#include <iostream>
#include <fstream>
#include <string>

inline bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}
