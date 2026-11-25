#include "alignment.hpp"
#include "multiple_alignment.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

int main() {
    const auto gapped = align_global("ACGT", "AGT");

    assert(gapped.first == "ACGT");
    assert(gapped.second == "A-GT");
    assert(gapped.score == 3);
    assert(gapped.matches == 3);
    assert(gapped.gaps == 1);
    assert(gapped.identity_percent() == 75.0);

    const auto mismatch = align_global("AC", "GT");

    assert(mismatch.score == -2);
    assert(mismatch.matches == 0);
    assert(mismatch.gaps == 0);
    assert(unique_pairs(0).empty());
    assert(unique_pairs(4) == (std::vector<std::pair<std::size_t, std::size_t>>{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}));

    const auto protein_pair = align_global("MA", "MTA");
    const auto codons = project_codons(protein_pair, "ATGGCTTAA", "ATGACTGCTTAA");

    assert(codons.first == "ATG---GCTTAA");
    assert(codons.second == "ATGACTGCTTAA");
    assert(codons.gaps == 3);

    const auto multiple = align_progressively({"AC", "ABC", "ADC"});
    assert(multiple.size() == 3);
    assert(multiple[0].size() == multiple[1].size());
    assert(multiple[1].size() == multiple[2].size());
    for (std::size_t i = 0; i < multiple.size(); ++i) {
        auto ungapped = multiple[i];
        ungapped.erase(std::remove(ungapped.begin(), ungapped.end(), '-'), ungapped.end());
        assert(ungapped == (std::vector<std::string>{"AC", "ABC", "ADC"})[i]);
    }
    const auto columns = summarize_columns(multiple);
    assert(columns.conserved == 2);
    assert(columns.with_gap == 1);

    const auto projected = project_multiple_codons(multiple,
        {"GCTTGTTAA", "GCTGGTTGTTAA", "GCTGATTGTTAA"});
    assert(projected[0].size() == 3 * (multiple[0].size() + 1));
    assert(projected[0].find("---") != std::string::npos);

    const auto synonymous = compare_codons(align_global("MS", "MS"),
                                           "ATGAGCTAA", "ATGAGTTAA");
    assert(synonymous.different_bases == 1);
    assert(synonymous.synonymous_codons == 1);
    assert(synonymous.nonsynonymous_codons == 0);
    assert(synonymous.gap_codons == 0);

    const std::vector<std::string> sequences{"ACGT", "AGT", "ACCT", "ACGA"};
    const auto serial = align_pairs(sequences, 1);
    const auto parallel = align_pairs(sequences, 4);

    assert(serial.size() == 6);

    for (std::size_t i = 0; i < serial.size(); ++i) {
        assert(serial[i].first_index == parallel[i].first_index);
        assert(serial[i].second_index == parallel[i].second_index);
        assert(serial[i].alignment.first == parallel[i].alignment.first);
        assert(serial[i].alignment.second == parallel[i].alignment.second);
        assert(serial[i].alignment.score == parallel[i].alignment.score);
    }

    std::cout << "alignment tests passed\n";
}
