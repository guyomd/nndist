#pragma once
#include <string>
#include <vector>
#include <map>
#include "nnanalysis.h"


bool readFromCSV(const std::string& filename, 
                 std::vector<std::vector<double>>& columns,
                 double minmag,
                 const std::vector<std::string> expected_headers,
                 char delimiter = ';');

bool writeToCSV(const std::string& filename, 
                const Hypocenters& events, 
                const std::vector<std::vector<double>>& values, 
                const std::vector<std::string> values_headers,
                char delimiter = ';');

std::map<std::string, std::string> readParametersFile(const std::string& filename);

