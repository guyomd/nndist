#pragma once
#include <vector>
#include <cstddef>
#include <valarray>

std::valarray<double> greatCircleDistance(double lat1, 
                                          double lon1, 
                                         std::valarray<double> lat2, 
                                         std::valarray<double> lon2);
std::valarray<double> horizontalDistance(double lat1, 
                                        double lon1,  
                                        std::valarray<double> lat2, 
                                        std::valarray<double> lon2);


// Core nearest neighbor computation (used internally)
std::vector<double> computeNearestNeighbor(
    const double lat,
    const double lon,
    const double dep,
    const double date,
    const std::valarray<double>& lats,
    const std::valarray<double>& lons,
    const std::valarray<double>& deps,
    const std::valarray<double>& dates,
    const std::valarray<double>& mags,
    char coords_unit,
    double w, double d, double p, double q
);

class Hypocenters {
public:
    void loadEvents(std::vector<std::vector<double>>& columns, 
                    size_t date_column_idx, 
                    char unit);
    void nearestNeighborDistance(double w, 
                                 double d, 
                                 double p, 
                                 double q);
    std::vector<std::vector<double>> decluster(double eta0, 
                                               double alpha0, 
                                               double w,
                                               double d,
                                               size_t npert, 
                                               double p = 0.5, 
                                               double q = 0.5);
    size_t size() const;

    // Data members (public for simplicity, but you may want to make them private with accessors)
    std::valarray<double> dates, lats, lons, mags, deps;
    std::valarray<double> nndist, R, T, parent_mag;
    std::vector<size_t> parent_idx;
    size_t nev = 0;
    char coords_unit = 'd';
};