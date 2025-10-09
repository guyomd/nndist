#include <vector>
#include <algorithm>  // for find_if()

void sortColumnsAsOne(std::vector<std::vector<double>>& columns, size_t col_index = 0) {
    if (columns.empty() || columns[0].empty()) return;
    size_t n = columns[col_index].size();
    size_t m = columns.size();

    // Create index vector [0, 1, ..., n-1]
    std::vector<size_t> idx(n);
    for (size_t i = 0; i < n; ++i) idx[i] = i;

    // Sort indices based on columns[0]
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        return columns[col_index][a] < columns[col_index][b];
    });

    // Apply permutation to all columns
    for (size_t col = 0; col < m; ++col) {
        std::vector<double> sorted_col(n);
        for (size_t row = 0; row < n; ++row) {
            sorted_col[row] = columns[col][idx[row]];
        }
        columns[col] = std::move(sorted_col);
    }
}


