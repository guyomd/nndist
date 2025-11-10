#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <random>
#include <valarray>
#include <omp.h>
#include <string>
#include "arrayutils.h"
#include "nnanalysis.h"
#include "progressbar.hpp"

double inf = std::numeric_limits<double>::infinity();

// Utility: Compute great-circle distance (in km)
std::valarray<double> greatCircleDistance(double lat1, double lon1, 
                                          std::valarray<double> lat2, 
                                          std::valarray<double> lon2) {
    constexpr double radius = 6371;
    double d2r = M_PI / 180.0;
    std::valarray<double> dlat = (lat2 - lat1) * d2r;
    std::valarray<double> dlon = (lon2 - lon1) * d2r;
    auto sin2term = [](std::valarray<double> dang) { return pow(sin(0.5 * dang), 2); };
    return 2.0 * radius * asin(sqrt(sin2term(dlat) + cos(lat1 * d2r) * cos(lat2 * d2r) * sin2term(dlon)));
}

// Utility: Compute horizontal distance (in km)
std::valarray<double>  horizontalDistance(double lat1, double lon1, 
                                          std::valarray<double> lat2, 
                                          std::valarray<double> lon2) {
    std::valarray<double>  dlat = lat2 - lat1;
    std::valarray<double>  dlon = lon2 - lon1;
    return sqrt(pow(dlon, 2) + pow(dlat, 2));
}


// Generic nearest neighbor computation (core logic)
std::vector<double> computeNearestNeighbor(const double lat,
                                           const double lon,
                                           const double dep,
                                           const double date,
                                           const std::valarray<double>& lats,
                                           const std::valarray<double>& lons,
                                           const std::valarray<double>& deps,
                                           const std::valarray<double>& dates,
                                           const std::valarray<double>& mags,
                                           char coords_unit,
                                           double w, double d, double p, double q) 
{
    size_t ne = lats.size();
    std::valarray<double> horiz, spatial, eta(inf, ne), R, T;
    std::vector<double> min_dist{inf, inf, inf, -1, -999.9}; 
    std::valarray<double> dt = date - dates;
    std::vector<size_t> ineg = indicesSE(dt, 0.0);

    if (ineg.size() == ne) {
        // When dates has no value smaller than date
        return min_dist; 
    } else {
        // Change negative dt to inf so as to ensure a maximum eta below, and to 
        // discard them from the minimum nearest-neighbor search:
        for (const size_t& i: ineg) {
            dt[i] = inf;
        }
    }
    horiz = (coords_unit == 'd')
        ? greatCircleDistance(lat, lon, lats, lons)
        : horizontalDistance(lat, lon, lats, lons);
    spatial = sqrt(pow(horiz, 2) + pow(deps - dep, 2));
    R = pow(spatial, d) * pow(10, -w * p * mags);
    T = dt * pow(10, -w * q * mags);
    eta = R * T;
    int imin = std::distance(std::begin(eta), std::min_element(std::begin(eta), std::end(eta)));
    min_dist[0] = eta[imin];
    min_dist[1] = R[imin];
    min_dist[2] = T[imin];
    min_dist[3] = imin + 1;  // indexing starts at 1
    min_dist[4] = mags[imin];
    return min_dist;
}



// Hypocenters methods

void Hypocenters::loadEvents(std::vector<std::vector<double>>& columns, 
                             size_t date_column_idx, char unit) {
    sortColumnsAsOne(columns, date_column_idx);
    nev = columns[0].size();
    if (columns.size() < 5) throw std::runtime_error("Not enough columns");
    if (not std::is_sorted(columns[0].begin(), columns[0].end())) {
        std::cerr << "input dates not in ascending order" << std::endl;
        exit(10); 
    }
    dates = std::valarray<double>(columns[0].data(), columns[0].size());
    lats  = std::valarray<double>(columns[1].data(), columns[1].size());
    lons  = std::valarray<double>(columns[2].data(), columns[2].size());
    mags  = std::valarray<double>(columns[3].data(), columns[3].size());
    deps = std::valarray<double>(columns[4].data(), columns[4].size());
    if (columns.size() > 5) {
        nndist = std::valarray<double>(columns[5].data(), columns[5].size());
        R = std::valarray<double>(columns[6].data(), columns[6].size());
        T = std::valarray<double>(columns[7].data(), columns[7].size());
        parent_idx.resize(nev);
        for (size_t k = 0; k < nev && k < columns[8].size(); ++k)
            parent_idx[k] = static_cast<size_t>(columns[8][k]);
        parent_mag = std::valarray<double>(columns[9].data(), columns[9].size());
    }
    coords_unit = unit;
}

void Hypocenters::nearestNeighborDistance(double w, double d, double p, double q) 
{
    nndist.resize(nev, inf);
    R.resize(nev, inf);
    T.resize(nev, inf);
    parent_idx.resize(nev, 1);
    parent_mag.resize(nev, -999.0);

    std::cout << ">> Compute nearest neighbor distances..." << std::endl;
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 1; i < nev; ++i) {
        auto min_dist = computeNearestNeighbor(lats[i], lons[i], deps[i], dates[i],
                                                lats[std::slice(0, i, 1)], 
                                                lons[std::slice(0, i, 1)], 
                                                deps[std::slice(0, i, 1)], 
                                                dates[std::slice(0, i, 1)], 
                                                mags[std::slice(0, i, 1)], 
                                                coords_unit, w, d, p, q);
        if (std::isinf(min_dist[0]) && (min_dist[4] == -999.0)) {
            min_dist[3] = i;  // set own event index
        }
        nndist[i] = min_dist[0];
        R[i] = min_dist[1];
        T[i] = min_dist[2];
        parent_idx[i] = min_dist[3];
        parent_mag[i] = min_dist[4];
    }
}

size_t Hypocenters::size() const {
    return nev;
}

std::vector<std::vector<double>> Hypocenters::decluster(double eta0, double alpha0, double w, 
                                                        double d, size_t npert, double p, double q,
                                                        std::string t_sampling_mode) 
{
    std::cout << ">> Applying nearest-neighbor declustering (Zaliapin & Ben-Zion, 2020)" << std::endl;
    std::cout << ">> Parameters: " << std:: endl;
    std::cout << "   -- w = " << w << std::endl;
    std::cout << "   -- d = " << d << std::endl;
    std::cout << "   -- p = " << p << std::endl;
    std::cout << "   -- q = " << q << std::endl;
    std::cout << "   -- eta0 = " << eta0 << std::endl;
    std::cout << "   -- alpha0 = " << alpha0 << std::endl;
    std::vector<size_t> i0 = indicesGE<double>(nndist, eta0);
    std::sort(i0.begin(), i0.end());
    size_t nm = i0.size(), idx;
    std::cout << "   -- " << nm << " events with eta0 >= " << eta0 << std::endl;
    double tmin = *std::min_element(std::begin(dates), std::end(dates));
    double tmax = *std::max_element(std::begin(dates), std::end(dates));
    std::valarray<double> m_k(nm), t_k(nm), x_k(nm), y_k(nm), z_k(nm);
    for (size_t i = 0; i < nm; ++i) {
        m_k[i] = mags[i0[i]];
        t_k[i] = dates[i0[i]];
        x_k[i] = lons[i0[i]];
        y_k[i] = lats[i0[i]];
        z_k[i] = deps[i0[i]];
    }
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::valarray<std::valarray<double>> nnreal(std::valarray<double>(inf, nev), npert);
    double dt;
    bool cond_x, cond_y;
    progressbar bar(npert);

    std::cout << "   -- " << npert; 
    std::cout << " random permutations of background events..." << std::endl;
    std::cout << "   -- randomization mode for origin times: " + t_sampling_mode << std::endl;
    #pragma omp parallel for schedule(dynamic)
    for (size_t k = 0; k < npert; ++k) {
        // Display progress:
        #pragma omp critical
            bar.update();
        
        // Each thread has its own random number generator
        thread_local std::mt19937 local_rng(std::random_device{}() + omp_get_thread_num());
        
        // Local thread variables
        std::valarray<double> local_m_k = m_k;  // local copy
        std::valarray<double> local_t_k = t_k;  // local copy
        

        // Random permutations of magnitudes (on local copy):
        std::shuffle(std::begin(local_m_k), std::end(local_m_k), local_rng); 
        
        if (t_sampling_mode == "synthetic") {
            // Draw uniform random occurrence times between tmin and tmax:
            for (auto& t : local_t_k)
                t = tmin + (tmax - tmin) * unif(local_rng);
        } 
        else if (t_sampling_mode == "permute") {
            // Random permutations of origin times (on local copy):
            std::shuffle(std::begin(local_t_k), std::end(local_t_k), local_rng);
        }
        
        for (size_t i = 1; i < nev; ++i) {
            std::valarray<bool> cond(nm);
            double dt;
            bool cond_x, cond_y;
            
            for (size_t j = 0; j < nm; ++j) {
                dt = dates[i] - local_t_k[j];
                cond_x = (lons[i] != x_k[j]);
                cond_y = (lats[i] != y_k[j]);
                cond[j] = (dt > 0.0) && (cond_x && cond_y);
            }
            
            std::vector<size_t> ivalid = indicesEQ<bool>(cond, true);
            size_t nval = ivalid.size();
            if (nval == 0) continue;
            std::valarray<double> x_l(nval);
            std::valarray<double> y_l(nval);
            std::valarray<double> z_l(nval);
            std::valarray<double> t_l(nval);
            std::valarray<double> m_l(nval); 
            for (size_t l = 0; l < nval; ++l) {
                idx = ivalid[l];
                x_l[l] = x_k[idx];
                y_l[l] = y_k[idx];
                z_l[l] = z_k[idx];
                t_l[l] = local_t_k[idx];
                m_l[l] = local_m_k[idx];
            }
            auto nnres = computeNearestNeighbor(lats[i], lons[i], deps[i], dates[i],
                                                y_l, x_l, z_l, t_l, m_l, 
                                                coords_unit, w, d, p, q);
            nnreal[k][i] = std::move(nnres[0]);  // NB: nnres = [nndist, R, T, parent_idx, parent_mag]
        }
    }
    std::cerr << std::endl;  // Add line break when progress bar is completed
    
    thread_local std::mt19937 rng(std::random_device{}());
    std::vector<std::vector<double>> results(nev, std::vector<double>(4, 0.0));
    results[0] = {0.0, 0, 0.0, inf};  // for first event
    double A0 = pow(10, alpha0);
    for (size_t i = 1; i < nev; ++i) {
        double count = 0.0, avg_lognnd = 0.0;
        for (size_t k = 0; k < npert; ++k) {
            if (std::isfinite(nnreal[k][i])) {
                count += 1.0;
                avg_lognnd += log10(nnreal[k][i]);
            }
        }
        avg_lognnd = (count > 0) ? avg_lognnd / count : inf;  // Average log10(nearest-neighbor distance)
        double prox = A0 * pow(10, log10(nndist[i]) - avg_lognnd);
        double prob = std::min(prox, 1.0);
        //double isbng = static_cast<double>(prob >= unif(rng));
        results[i] = {prob, prox, pow(10, avg_lognnd)};  // [prob_bgnd, norm. prox., avg_nnd]
    }
    return results;
}


std::vector<double> Hypocenters::extractBgndEventTimes(const std::vector<std::vector<double>>& decluster_results) 
{
     /*
     * Extract background event times by random thinning from original catalogue and declustering results
     * 
     * @param dates All event occurrence times
     * @param decluster_results Results from decluster method (each element: [prob_bgnd, norm_prox, avg_nnd])
     * @return Vector of occurrence times for events where is_bgnd == 1 (True)
     */
    std::vector<double> background_times;
    size_t n = std::min(dates.size(), decluster_results.size());
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    thread_local std::mt19937 rng(std::random_device{}());
    
    for (size_t i = 0; i < n; ++i) {
        // Random thinning: Keep event if prob_bgnd[i] >= rand. sample U(0, 1):
        if (decluster_results[i][0] >= unif(rng)) {
            background_times.push_back(dates[i]);
        }
    }
    
    // Sort the times
    std::sort(background_times.begin(), background_times.end());
    return background_times;
}


StatTestResult Hypocenters::testBgndStationarityBrownZhao(const std::vector<double>& background_times,
                                                          int k,
                                                          double alpha) 
{
    if (background_times.size() < 2) {
        std::cout << "Warning: Insufficient background events (" << background_times.size() 
                  << ") for Brown-Zhao test." << std::endl;
        StatTestResult result;
        result.test_name = "Brown-Zhao (2002) test";
        result.statistic = 0.0;
        result.p_value = 1.0;
        result.is_stationary = true;
        return result;
    }
    
    return PoissonStationarityTests::brownZhaoTest(background_times, k, alpha);
}


StatTestResult Hypocenters::testBgndStationarityKS(const std::vector<double>& background_times,
                                                   double alpha) 
{
    if (background_times.size() < 2) {
        std::cout << "Warning: Insufficient background events (" << background_times.size() 
                  << ") for Kolmogorov-Smirnov test." << std::endl;
        StatTestResult result;
        result.test_name = "Kolmogorov-Smirnov test";
        result.statistic = 0.0;
        result.p_value = 1.0;
        result.is_stationary = true;
        return result;
    }
    
    return PoissonStationarityTests::kolmogorovSmirnovTest(background_times, alpha);
}


void Hypocenters::performStationarityTests(const std::vector<std::vector<double>>& decluster_results,
                                           double alpha, int n) 
{
    std::vector<double> bz1(n), bz2(n), ks(n);
    StatTestResult bz1_result, bz2_result, ks_result;

    std::cout << "\n### STATIONARITY TESTS FOR BACKGROUND EVENTS ###" << std::endl;
    std::cout << "-- Based on " << n << " random background catalogues" << std::endl;
    
    for (size_t i = 0; i < n; ++i) {
	// Extract a realization of background events:
	std::vector<double> background_times = extractBgndEventTimes(decluster_results);
	
	/*
        std::cout << "Total events: " << dates.size() << std::endl;
	std::cout << "Background events: " << background_times.size() << std::endl;
	*/
        if (background_times.size() > 0) {
	    /*
            std::cout << "Background event rate: " << 
			 static_cast<double>(background_times.size()) / dates.size() * 100.0 
		      << "%" << std::endl;
	    std::cout << "Time span: " << background_times.front() << " to " 
		      << background_times.back() << std::endl;
            */
	}
	
	if (background_times.size() < 2) {
	    std::cout << "\nInsufficient background events for statistical tests." << std::endl;
	    return;
	}
	
	// Perform Brown-Zhao tests:
	bz1_result = testBgndStationarityBrownZhao(background_times, 10, alpha);
        bz1[i] = bz1_result.p_value;
	bz2_result = testBgndStationarityBrownZhao(background_times, 100, alpha);
        bz2[i] = bz2_result.p_value;

	// // Perform Kolmogorov-Smirnov test
	ks_result = testBgndStationarityKS(background_times, alpha);
        ks[i] = ks_result.p_value;
    }
    
    if (n > 1) {
	// Compute median, 2.5% and 97.5% percentiles for each test:
	bz1_result.p_range[0] = getPercentile(bz1, 2.5);
	bz1_result.p_range[1] = getPercentile(bz1, 50.0);
	bz1_result.p_range[2] = getPercentile(bz1, 97.5);
        bz1_result.p_value = bz1_result.p_range[1]; // Replace by median p-value
        bz1_result.is_stationary = (bz1_result.p_value > alpha);  // Update result
	bz1_result.printTestResults();

	bz2_result.p_range[0] = getPercentile(bz2, 2.5);
	bz2_result.p_range[1] = getPercentile(bz2, 50.0);
	bz2_result.p_range[2] = getPercentile(bz2, 97.5);
        bz2_result.p_value = bz2_result.p_range[1]; 
        bz2_result.is_stationary = (bz2_result.p_value > alpha); 
	bz2_result.printTestResults();

	ks_result.p_range[0] = getPercentile(ks, 2.5);
	ks_result.p_range[1] = getPercentile(ks, 50.0);
	ks_result.p_range[2] = getPercentile(ks, 97.5);
        ks_result.p_value = ks_result.p_range[1]; 
        ks_result.is_stationary = (ks_result.p_value > alpha); 
	ks_result.printTestResults();
    } 
    else {
	bz1_result.printTestResults();
	bz2_result.printTestResults();
	ks_result.printTestResults();
    }
    
	
}
