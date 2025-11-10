#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <cstdlib>
#include <iterator>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/kolmogorov_smirnov.hpp>
#include "nnstats.h"



double getPercentile(const std::vector<double>& samples, double percentile)
{
    std::vector<double> values(samples);  

    // Define index of percentile position based on vector size:
    const std::size_t pos = (percentile / 100) * values.size();   
    auto pth = values.begin() + pos;
    std::nth_element(values.begin(), pth, values.end());  
    return values[pos];

}


double chiSquarePValue(double testStat, int degreesOfFreedom)
{
    boost::math::chi_squared_distribution<double> chiSquared(degreesOfFreedom);
    double cdfValue = boost::math::cdf(chiSquared, testStat);
    double pValue = 1.0 - cdfValue;
    return pValue;
}


double ksPValue(double testStat, double n)
{
    boost::math::kolmogorov_smirnov_distribution<double> ks_dist(n);
    double cdfValue = boost::math::cdf(ks_dist, testStat);
    double pValue = 1.0 - cdfValue;
    return pValue;
}


StatTestResult PoissonStationarityTests::brownZhaoTest(const std::vector<double>& occurrence_times, 
                                                       int k,
                                                       double alpha) 
{
    StatTestResult result;
    result.test_name = "Brown-Zhao (2002) test: k=" + std::to_string(k) ;

    std::vector<double> yk(k);
    double ymean = 0;
    double bzstat = 0;

    // Clean and sort the data
    std::vector<double> times = sortAndClean(occurrence_times);
    int n = times.size();
    double tmin = *std::min_element(times.begin(), times.end());
    double tmax = *std::max_element(times.begin(), times.end());
    double tsub = (tmax - tmin) / k;
   
    // Compute numbers in each (equal-duration sub-interval):
    for (size_t i = 0; i < k; ++i) {
        double tb = tmin + i * tsub;
        double te = tmin + (i + 1) * tsub;
        int count = 0;
        auto lower = std::lower_bound(times.begin(), times.end(), tb);
        auto upper = std::upper_bound(times.begin(), times.end(), te);
        for (auto it = lower; it != upper; it++) {
            count += 1;
        }
        yk[i] = std::sqrt(count + 3 / 8);
        ymean += yk[i];
    }
    ymean /= k;

    // Compute Brown & Zhao (2002) test statistics:
    for (size_t i = 0; i < k; ++i) {
        bzstat += 4 * std::pow(yk[i] - ymean, 2);
    }
    result.statistic = bzstat;
    result.p_value = chiSquarePValue(bzstat, k - 1);
    result.is_stationary = (result.p_value > alpha);
    return result;
}


StatTestResult PoissonStationarityTests::kolmogorovSmirnovTest(const std::vector<double>& occurrence_times,
                                                               double alpha) 
{
    // Clean and sort the data
    std::vector<double> times = sortAndClean(occurrence_times);
    int n = times.size();
    double tmin = *std::min_element(times.begin(), times.end());
    double tmax = *std::max_element(times.begin(), times.end());

    StatTestResult result;
    result.test_name = "One-sided Kolmogorov-Smirnov test (" + std::to_string(n) + " samples)";
    
    // Compute transformed times:
    std::vector<double> u(n), unif(n);  
    for (size_t i = 0; i < n;  ++i) {
        u[i] = (times[i] - tmin) / (tmax - tmin);  
    }
    
    // Compute empirical and theoretical (uniform) CDF:
    std::vector<double> u_CDF(n), theo_CDF(n), diff(n);
    double nd = static_cast<double>(n);
    for (size_t i = 0; i < n;  ++i) {
        u_CDF[i] = static_cast<double>(i) / (nd - 1.0);
        theo_CDF[i] = u[i];
        diff[i] = std::abs(u_CDF[i] - theo_CDF[i]);
    }
 
    result.statistic = *std::max_element(diff.begin(), diff.end());
    result.p_value = ksPValue(result.statistic, static_cast<double>(n));
    result.is_stationary = (result.p_value > alpha);
    return result;
}


std::vector<double> PoissonStationarityTests::sortAndClean(const std::vector<double>& data) {
    std::vector<double> cleaned = data;
    std::sort(cleaned.begin(), cleaned.end());
    
    // Remove duplicates
    auto it = std::unique(cleaned.begin(), cleaned.end());
    cleaned.resize(std::distance(cleaned.begin(), it));
    
    return cleaned;
}


void StatTestResult::printTestResults() 
{
    std::cout << "\n=== " << test_name << " ===" << std::endl;
    if ((p_range[0] == -1.0) and (p_range[1] == -1.0) and (p_range[2] == -1.0)) {
        std::cout << "Test statistic: " << statistic << std::endl;
        std::cout << "P-value: " << p_value << std::endl;
    }
    else {
        std::cout << "P-value:  median = " << p_range[1] << ", range = [" << p_range[0] 
                  << ", " << p_range[2] << "]  (Q2.5%, Q97.5%)" << std::endl;
    }
    std::cout << "Result: ";
    if (is_stationary) {
        std::cout << "FAIL TO REJECT null hypothesis (process appears stationary)" << std::endl;
    } else {
        std::cout << "REJECT null hypothesis (process appears non-stationary)" << std::endl;
    }
    std::cout << std::endl;
}
