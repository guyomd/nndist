#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <ctime>  // for std::clock_t and std::clock()
#include <iomanip>
#include <stdexcept>
#include <omp.h>  // parallelization
#include "strutils.h"
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
    std::cout << ">> Using configuration file: " << cfgfile << std::endl;
    std::map<std::string, std::string> params = readParametersFile(cfgfile);
    std::string infile = params["output_file"];
    std::string prefix = infile.substr(0, infile.find('.'));
    std::string suffix = infile.substr(infile.find('.'));
    std::string outfile = prefix + "_declust" + suffix;
    
    // Not sure we need all of these:
    char unit = params["unit_for_geog_coordinates_[km,degrees]"][0];
    double d = std::stod(params["fractal_dimension"]);
    double p = std::stod(params["parameter_p"]);
    double q = std::stod(params["parameter_q"]);
    double w = std::stod(params["parameter_w"]);
    double eta0 = std::stod(params["parameter_eta0"]);
    std::vector<double> alpha0_values;
    try {
        alpha0_values = parseAlpha0Values(params["parameter_alpha0"]);
    } catch (const std::exception& ex) {
        std::cerr << "Error parsing parameter_alpha0: " << ex.what() << std::endl;
        return 1;
    }
    if (alpha0_values.empty()) {
        std::cerr << "Error: parameter_alpha0 must specify at least one value." << std::endl;
        return 1;
    }
    size_t npert = std::stoul(params["parameter_npert"]);
    std::string t_mode = params["time_perturbation_mode_[synthetic,permute]"];
    int ntests = std::stoi(params["nb_realizations_for_stats_tests"]);

    std::clock_t start_time = std::clock();
    std::vector<std::vector<double>> columns;
    size_t date_col_idx = 0;
    Hypocenters events;
    std::vector<std::vector<double>> results;
    std::vector<std::string> expected_headers = {"floating_date", "latitude", "longitude", "magnitude", "depth", 
                                                 "nn_distance", "R", "T", "index_ancestor", "mag_ancestor"}; 
    if (readFromCSV(infile, columns, (double) -999.99, expected_headers, ';')) {
        events.loadEvents(columns, date_col_idx, unit);
        std::vector<std::string> results_headers = {"p_bgnd", "norm_prox", "avg_nn_distance"};
        auto all_results = events.decluster(eta0, alpha0_values, w, d, npert, p, q, t_mode);
        for (size_t ai = 0; ai < alpha0_values.size(); ++ai) {
            double alpha0 = alpha0_values[ai];
            const auto& results = all_results[ai];
            std::cout << "\n### NEAREST-NEIGHBOR DECLUSTERING COMPUTED FOR ALPHA0 = " << alpha0 << std::endl;

            if (ntests > 0) {
                double test_alpha = 0.05;
                events.performStationarityTests(results, test_alpha, ntests);
            }

            std::string alpha0_outfile = outfile;
            if (alpha0_values.size() > 1) {
                alpha0_outfile = prefix + "_declust_alpha0_" + formatAlpha0(alpha0) + suffix;
            }

            if (writeToCSV(alpha0_outfile, events, results, results_headers, ';')) {
                std::cout << ">> Results saved in file " << alpha0_outfile << std::endl;
            } else {
                std::cerr << "Error writing results to file " << alpha0_outfile << std::endl;
            }
        }
        std::clock_t end_time = std::clock();
        double elapsed_seconds = static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cout << ">> Total elapsed CPU time: " << std::fixed << std::setprecision(2) << elapsed_seconds << " s" << std::endl;
    }
    return 0;
}
