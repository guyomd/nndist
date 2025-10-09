#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
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
 
    // Print first line using mixed fixed-width and scientific notations:
    outfile << std::fixed 
            << std::setprecision(12) << events.dates[0] << sep
            << std::setprecision(6) << events.lats[0] << sep
            << std::setprecision(6) << events.lons[0] << sep
            << std::setprecision(2) << events.mags[0] << sep
            << std::setprecision(3) << events.deps[0] << sep
            << std::scientific
            << std::setprecision(2) << events.nndist[0] << sep
            << std::setprecision(2) << events.R[0] << sep
            << std::setprecision(2) << events.T[0] << sep
            << std::fixed
            << std::setprecision(0) << events.parent_idx[0] << sep
            << std::setprecision(2) << events.parent_mag[0];
    if (nval > 1) {
        outfile << sep 
                << std::fixed 
                << std::setprecision(2) << values[0][0] << sep  // 1.O
                << std::setprecision(0) << values[0][1] << sep  // 1
                << std::scientific
                << std::setprecision(2) << values[0][2] << sep  // 1.0
                << std::setprecision(2) << values[0][3];        // 1E99        
    }
    outfile << std::endl;
    
    // Then, print other lines using fixed-width notation:
    outfile << std::fixed;  
    for (size_t i = 1; i < events.nev; ++i) {
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
                    << std::setprecision(0) << values[i][1] << sep  // Background event flag (0/1)
                    << std::setprecision(14) << values[i][2] << sep  // Normalized proximity 
                    << std::setprecision(14) << values[i][3];        // Avg. nearest-neighbor distance over permutations
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