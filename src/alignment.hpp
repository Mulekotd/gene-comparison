#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

struct Alignment {
    std::string first;
    std::string second;
    std::size_t matches = 0;
    std::size_t gaps = 0;
    int score = 0;

    double identity_percent() const;
};

struct PairResult {
    std::size_t first_index;
    std::size_t second_index;
    Alignment alignment;
};

Alignment align_global(const std::string& first, const std::string& second, int match = 2, int mismatch = -1, int gap = -3);
Alignment project_codons(const Alignment& protein_alignment, const std::string& first_cds, const std::string& second_cds);
std::vector<std::pair<std::size_t, std::size_t>> unique_pairs(std::size_t sequence_count);
std::vector<PairResult> align_pairs(const std::vector<std::string>& sequences, std::size_t thread_count);
