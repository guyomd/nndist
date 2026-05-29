#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "strutils.h"
#include "nnanalysis.h"


// Reads columns from a CSV input file.
// Each column is loaded into a vector<double> in the output vector.
bool readFromCSV(const std::string& filename, 
                 std::vector<std::vector<double>>& columns, 
                 double minmag, 
                 const std::vector<std::string> expected_headers,
                 char delimiter = ';') {
    
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "[readFromCSV] Error opening file: " << filename << std::endl;
        return false;
    }
    std::cout << ">> Loading data from file: " << filename << std::endl;

    std::string line;
    size_t num_columns = 0;
    int line_num = 0;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string value;
        //static const std::vector<std::string> expected_headers = {"floating_date", "latitude", "longitude", "magnitude", "depth"};
        static bool header_processed = false;
        std::vector<double> row_values;
        
        line_num += 1;

        // Header processing block
        if (!header_processed) {
            std::istringstream header_stream(line);
            std::string header;
            std::vector<std::string> headers;
            while (std::getline(header_stream, header, delimiter)) {
                header = trim(header);
                if (header.size() > 0) {
                    headers.push_back(trim(header));
                }
            }
       
            if (!std::equal(headers.begin(), headers.end(), expected_headers.begin(), expected_headers.end())) {
                std::cerr << "[readFromCSV] Header mismatch. Expected columns: ";
                for (size_t c = 0; c < expected_headers.size(); c++) {
                    std::cerr << "[" << expected_headers[c] << "] ";
                }
                std::cerr << std::endl;
                return false;
            }
            header_processed = true;
            continue; // Skip header line
        }

        // Row processing block
        while (std::getline(iss, value, delimiter)) {
            try {
                value = trim(value);
                if (value.size() > 0) {
                    row_values.push_back(std::stod(value));
                }  
            } catch (...) {
                std::cerr << "[readFromCSV] l. " << line_num << ": Unexpected format for value " << value << std::endl;
                return false;
            }
        }
        if (columns.empty()) {
            num_columns = row_values.size();
            columns.resize(num_columns);
        }
        if (row_values[3] >= minmag) {  // store only if magnitude is above minimum threshold
            for (size_t i = 0; i < num_columns && i < row_values.size(); i++) {
                columns[i].push_back(row_values[i]);
            }
        }
    }
    infile.close();
    return true;
}



bool writeToCSV(const std::string& filename, const Hypocenters& events, 
                const std::vector<std::vector<double>>& values, 
                const std::vector<std::string> values_headers,
                char delimiter) 
{
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "[writeToCSV] Error opening file for writing: " << filename << std::endl;
        return false;
    }

    // Write header
    std::string sep = std::string(1, delimiter) + " ";
    size_t nval = values_headers.size();
    outfile << "floating_date" << sep
            << "latitude" << sep
            << "longitude" << sep
            << "magnitude" << sep
            << "depth" << sep
            << "nn_distance" << sep
            << "R" << sep
            << "T" << sep
            << "index_ancestor" << sep
            << "mag_ancestor";
    for (size_t i = 0; i < nval; ++i) {
        outfile << sep << values_headers[i];
    }
    outfile << std::endl;

    if ((nval > 0) && (values.size() != events.nev)) {
        std::cerr << "[writeToCSV] Mismatch between input values and events count: " << values[0].size() << " vs. " << events.nev << std::endl;
        return false;
    }
 
    // Then, print other lines using fixed-width notation:
    outfile << std::fixed;  
    for (size_t i = 0; i < events.nev; ++i) {
        outfile << std::setprecision(12) << events.dates[i] << sep
                << std::setprecision(6) << events.lats[i] << sep
                << std::setprecision(6) << events.lons[i] << sep
                << std::setprecision(2) << events.mags[i] << sep
                << std::setprecision(3) << events.deps[i] << sep
                << std::setprecision(14) << events.nndist[i] << sep
                << std::setprecision(14) << events.R[i] << sep
                << std::setprecision(14) << events.T[i] << sep
                << std::setprecision(0) << events.parent_idx[i] << sep
                << std::setprecision(2) << events.parent_mag[i];
        if (nval > 1) {
            outfile << sep 
                    << std::setprecision(6) << values[i][0] << sep  // prob. to be a background event
                    << std::setprecision(14) << values[i][1] << sep  // Normalized proximity 
                    << std::setprecision(14) << values[i][2];        // Avg. nearest-neighbor distance over permutations
        }
        outfile << std::endl;
    }
    outfile.close();
    return true;
}


// Read parameters file and return a dictionary of paremeter values (key-value pairs)
std::map<std::string, std::string> readParametersFile(const std::string& filename) {
    std::map<std::string, std::string> params;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error opening parameters file: " << filename << std::endl;
        return params;
    }
    std::string line;
    while (std::getline(infile, line)) {
        // Ignore empty lines and comments starting wth "#" symbol
        if (line.empty() || line[0] == '#') continue;
        auto pos = line.find(':');
        if (pos == std::string::npos) continue; // Invalid line
        std::string key = trim(line.substr(0, pos));
        std::string value = trim(line.substr(pos + 1));
        if (!key.empty()) {
            params[key] = value;
        }
    }
    infile.close();
    return params;
}

std::vector<double> parseAlpha0Values(const std::string& alpha0_param) {
    std::vector<double> values;
    std::string trimmed = trim(alpha0_param);
    if (trimmed.empty()) {
        return values;
    }

    if (trimmed.find(',') != std::string::npos) {
        auto tokens = splitString(trimmed, ',');
        for (const auto& token : tokens) {
            values.push_back(std::stod(token));
        }
        return values;
    }

    size_t colon_count = std::count(trimmed.begin(), trimmed.end(), ':');
    if (colon_count == 2) {
        auto parts = splitString(trimmed, ':');
        if (parts.size() != 3) {
            throw std::runtime_error("Invalid alpha0 range syntax: " + alpha0_param);
        }
        double start = std::stod(parts[0]);
        double step = std::stod(parts[1]);
        double end = std::stod(parts[2]);
        if (step == 0.0) {
            throw std::runtime_error("Alpha0 range step must not be zero.");
        }
        if ((step > 0.0 && start > end) || (step < 0.0 && start < end)) {
            throw std::runtime_error("Alpha0 range bounds are inconsistent with step direction.");
        }
        double value = start;
        const double eps = std::abs(step) * 1e-9;
        while ((step > 0.0 && value <= end + eps) || (step < 0.0 && value >= end - eps)) {
            values.push_back(value);
            value += step;
        }
        if (!values.empty() && std::abs(values.back() - end) > eps) {
            values.push_back(end);
        }
        return values;
    }

    values.push_back(std::stod(trimmed));
    return values;
}
