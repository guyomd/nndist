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
    std::string infile = params["output_file"];
    std::string prefix = infile.substr(0, infile.find('.'));
    std::string suffix = infile.substr(infile.find('.'));
    std::string outfile = prefix + "_declust" + suffix;
    
    // Not sure we need all of these:
    char unit = params["unit_for_geog_coordinates_[km,degrees]"][0];
    double d = std::stod( params["fractal_dimension"]);
    double p = std::stod( params["parameter_p"]);
    double q = std::stod( params["parameter_q"]);
    double w = std::stod( params["parameter_w"]);
    double eta0 = std::stod( params["parameter_eta0"]);
    double alpha0 = std::stod( params["parameter_alpha0"]);
    size_t npert = std::stod( params["parameter_npert"]);
    std::string t_mode = params["time_perturbation_mode_[synthetic,permute]"];
    int ntests = std::stoi(params["nb_realizations_for_stats_tests"]);

    // Load catalogue and nearest-neighbor distances:
    std::vector<std::vector<double>> columns;
    size_t date_col_idx = 0;
    Hypocenters events;
    std::vector<std::vector<double>> results;
    std::vector<std::string> expected_headers = {"floating_date", "latitude", "longitude", "magnitude", "depth", 
                                                 "nn_distance", "R", "T", "index_ancestor", "mag_ancestor"}; 
    if (readFromCSV(infile, columns, (double) -999.99, expected_headers, ';')) {
        events.loadEvents(columns, date_col_idx, unit);
        std::time_t c_start = std::time(nullptr); // NB: time resolution: 1 s.
        results = events.decluster(eta0, alpha0, w, d, npert, p, q, t_mode);  // default: p = q = 0.5
        std::time_t c_end = std::time(nullptr);
        double time_elapsed_s = std::difftime(c_end, c_start);
        std::cout << "CPU time used ~ " << time_elapsed_s << " s." << std::endl;
        
        // Perform stationarity tests on background events:
        double test_alpha = 0.05;
        events.performStationarityTests(results, test_alpha, ntests);
    }
    
    // Write results to output file:
    std::vector<std::string> results_headers = {"p_bgnd", "norm_prox", "avg_nn_distance"};
    if (writeToCSV(outfile, events, results, results_headers, ';')) {
        std::cout << ">> Results saved in file " << outfile << std::endl;
    }
    return 0;
}
