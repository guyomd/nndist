#pragma once
#include <vector>
#include <valarray>
#include <algorithm>
#include <iterator>

void sortColumnsAsOne(std::vector<std::vector<double>>& columns, size_t col_index = 0);

template <typename T>
std::vector<size_t> indicesGE(std::vector<T> v, double value)
// Returns a vector of indices matching the input logical condition
{
    std::vector<size_t> results;
    results.reserve(v.size());
    auto it = std::find_if(std::begin(v), std::end(v), [value](const T& elem){ return elem >= value; });
    while (it != std::end(v)) {
        results.emplace_back(std::distance(std::begin(v), it));
        it = std::find_if(std::next(it), std::end(v), [value](const T& elem){ return elem >= value; });
    }
    return results;
}

template <typename T>
std::vector<size_t> indicesGE(std::valarray<T> v, double value)
// Returns a vector of indices matching the input logical condition
{
    std::vector<size_t> results;
    results.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] >= value) {
            results.push_back(i);
        }
    }
    return results;
}

template <typename T>
std::vector<size_t> indicesGT(std::vector<T> v, double value)
// Returns a vector of indices matching the input logical condition
{
    std::vector<size_t> results;
    results.reserve(v.size());
    auto it = std::find_if(std::begin(v), std::end(v), [value](const T& elem){ return elem > value; });
    while (it != std::end(v)) {
        results.emplace_back(std::distance(std::begin(v), it));
        it = std::find_if(std::next(it), std::end(v), [value](const T& elem){ return elem > value; });
    }
    return results;
}

template <typename T>
std::vector<size_t> indicesGT(std::valarray<T> v, double value)
// Returns a vector of indices matching the input logical condition
{
    std::vector<size_t> results;
    results.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] > value) {
            results.push_back(i);
        }
    }
    return results;
}

template <typename T>
std::vector<size_t> indicesSE(std::vector<T> v, double value)
// Returns a vector of indices matching the input logical condition
{
    std::vector<size_t> results;
    results.reserve(v.size());
    auto it = std::find_if(std::begin(v), std::end(v), [value](const T& elem){ return elem <= value; });
    while (it != std::end(v)) {
        results.emplace_back(std::distance(std::begin(v), it));
        it = std::find_if(std::next(it), std::end(v), [value](const T& elem){ return elem <= value; });
    }
    return results;
}

template <typename T>
std::vector<size_t> indicesSE(std::valarray<T> v, double value)
// Returns a vector of indices matching the input logical condition
{
    std::vector<size_t> results;
    results.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] <= value) {
            results.push_back(i);
        }
    }
    return results;
}

template <typename T>
std::valarray<T> valuesSE(const std::valarray<T>& v, double value)
// Returns a valarray of values smaller than or equal to the given value
{
    std::vector<T> results;
    results.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] <= value) {
            results.push_back(v[i]);
        }
    }
    std::valarray<T> output(results.data(), results.size());
    return output;
}

template <typename T>
std::vector<size_t> indicesEQ(std::vector<bool> v, T value)
// Returns a vector of indices matching the input value
{
    std::vector<size_t> results;
    results.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] == value) {
            results.push_back(i);
        }
    }
    return results;
}

template <typename T>
std::vector<size_t> indicesEQ(std::valarray<bool> v, T value)
// Returns a vector of indices matching the input value
{
    std::vector<size_t> results;
    results.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] == value) {
            results.push_back(i);
        }
    }
    return results;
}