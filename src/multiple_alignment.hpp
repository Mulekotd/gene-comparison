#pragma once

#include "alignment.hpp"

#include <cstddef>
#include <string>
#include <vector>

struct ColumnSummary {
    std::size_t conserved = 0;
    std::size_t substituted = 0;
    std::size_t with_gap = 0;
};

struct CodonEffects {
    std::size_t different_bases = 0;
    std::size_t synonymous_codons = 0;
    std::size_t nonsynonymous_codons = 0;
    std::size_t gap_codons = 0;
    bool different_stop = false;
};

std::vector<std::string> align_progressively(const std::vector<std::string>& proteins);
std::vector<std::string> project_multiple_codons(const std::vector<std::string>& aligned_proteins, const std::vector<std::string>& coding_sequences);
ColumnSummary summarize_columns(const std::vector<std::string>& aligned_sequences);
CodonEffects compare_codons(const Alignment& aligned_proteins, const std::string& first_cds, const std::string& second_cds);
