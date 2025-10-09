#include <iostream>  
#include <vector>
#include <string>
#include <map>
#include <ctime>  // for std::clock_t and std::clock()
#include <omp.h>  // parallelization
#include "fileutils.h"
#include "nnanalysis.h"



// Returns the filename from command line arguments or a default if not provided
std::string getFilenameFromArgs(int argc, char* argv[], const std::string& defaultFilename = "config.txt") {
    if (argc > 1) {
        return argv[1];
    }
    return defaultFilename;
}



int main(int argc, char* argv[]) {

    // Parse parameters:
    std::string cfgfile = getFilenameFromArgs(argc, argv);
    std::cout << ">> Using configuration file: " << cfgfile << std::endl; // Example call and output
    std::map<std::string, std::string> params = readParametersFile(cfgfile);
    std::string catalog = params["input_catalogue"];
    std::string outfile = params["output_file"];
    double minmag = std::stod( params["minimum_magnitude"] );
    char unit = params["unit_for_geog_coordinates_[km,degrees]"][0];
    double d = std::stod( params["fractal_dimension"]);
    double p = std::stod( params["parameter_p"]);
    double q = std::stod( params["parameter_q"]);
    double w = std::stod( params["parameter_w"]);

    std::vector<std::vector<double>> columns;
    size_t date_col_idx = 0;
    Hypocenters events;
    std::vector<std::vector<double>> nn_values;
    
    struct timespec start, finish;
    double elapsed;

    // Load catalogue and compute nearest-neighbor distance:
    std::vector<std::string> expected_headers = {"floating_date", "latitude", "longitude", "magnitude", "depth"};
    if (readFromCSV(catalog, columns, minmag, expected_headers, ';')) {
        events.loadEvents(columns, date_col_idx, unit);
        std::cout << ">> Found " << events.size() << " events in catalogue with M >= " << minmag << std::endl;
        std::time_t c_start = std::time(nullptr); // NB: time resolution: 1 s.
        events.nearestNeighborDistance(w, d, p, q);
        std::time_t c_end = std::time(nullptr);
        double time_elapsed_s = std::difftime(c_end, c_start);
        std::cout << "CPU time used ~ " << time_elapsed_s << " s." << std::endl;
    }
    
    // Write results to output file:
    std::vector<std::vector<double>> fake_values(0, std::vector<double>(0));  //unused
    std::vector<std::string> fake_headers = {};  // when empty, function writeToCSV only writes Hypocenters object attributes
    if (writeToCSV(outfile, events, fake_values, fake_headers, ';')) {
        std::cout << ">> Results saved in file " << outfile << std::endl;
    }
    return 0;
}
