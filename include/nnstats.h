#include <vector>
#include <string>


class StatTestResult {
    public:
        double statistic;    // Test statistic value
        double p_value;      // P-value of the test
        bool is_stationary;  // True if null hypothesis (stationarity) is not rejected
        std::string test_name; // Name of the test performed
        void printTestResults();
};

class PoissonStationarityTests {
public:
    static StatTestResult brownZhaoTest(const std::vector<double>& occurrence_times, 
                                        int k,
                                        double alpha = 0.05);
    static StatTestResult kolmogorovSmirnovTest(const std::vector<double>& occurrence_times,
                                                double alpha = 0.05);
    static std::vector<double> sortAndClean(const std::vector<double>& data);

};

double chiSquarePValue(double testStat, int degreesOfFreedom);
double ksPValue(double testStat, double n);

