#pragma once
#include <string>
#include <vector>

std::string trim(const std::string& str, 
                 const std::string& whitespace = " \t");

std::vector<std::string> splitString(const std::string& text, char delimiter);

std::string formatAlpha0(double alpha0);
