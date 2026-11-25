#include "utils.hpp"

bool same_results(const std::vector<PairResult>& first, const std::vector<PairResult>& second) {
    if (first.size() != second.size())
        return false;

    for (std::size_t i = 0; i < first.size(); ++i) {
        if (first[i].first_index != second[i].first_index ||
            first[i].second_index != second[i].second_index ||
            first[i].alignment.score != second[i].alignment.score ||
            first[i].alignment.first != second[i].alignment.first ||
            first[i].alignment.second != second[i].alignment.second) return false;
    }

    return true;
}
