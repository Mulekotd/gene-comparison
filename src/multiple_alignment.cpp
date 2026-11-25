#include "multiple_alignment.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>

namespace {

constexpr int match_score = 2;
constexpr int mismatch_score = -1;
constexpr int gap_score = -3;

int pair_score(char first, char second) {
    if (first == '-' && second == '-')
        return 0;

    if (first == '-' || second == '-')
        return gap_score;

    return first == second ? match_score : mismatch_score;
}

std::size_t residue_count(const std::string& aligned) {
    return static_cast<std::size_t>(std::count_if(aligned.begin(), aligned.end(),
                                                   [](char ch) { return ch != '-'; }));
}

std::vector<std::string> add_sequence(const std::vector<std::string>& profile, const std::string& sequence) {
    const std::size_t profile_rows = profile.size();
    const std::size_t rows = profile.front().size() + 1;
    const std::size_t columns = sequence.size() + 1;

    std::vector<std::uint8_t> direction(rows * columns, 0);
    std::vector<int> previous(columns), current(columns);

    for (std::size_t j = 1; j < columns; ++j) {
        previous[j] = previous[j - 1] + gap_score * static_cast<int>(profile_rows);
        direction[j] = 2;
    }

    for (std::size_t i = 1; i < rows; ++i) {
        int gap_against_profile = 0;

        for (const auto& existing : profile) {
            gap_against_profile += pair_score(existing[i - 1], '-');
        }

        current[0] = previous[0] + gap_against_profile;
        direction[i * columns] = 1;

        for (std::size_t j = 1; j < columns; ++j) {
            int residue_against_profile = 0;

            for (const auto& existing : profile) {
                residue_against_profile += pair_score(existing[i - 1], sequence[j - 1]);
            }

            const int diagonal = previous[j - 1] + residue_against_profile;
            const int up = previous[j] + gap_against_profile;
            const int left = current[j - 1] + gap_score * static_cast<int>(profile_rows);

            if (diagonal >= up && diagonal >= left) {
                current[j] = diagonal;
                direction[i * columns + j] = 0;
            } else if (up >= left) {
                current[j] = up;
                direction[i * columns + j] = 1;
            } else {
                current[j] = left;
                direction[i * columns + j] = 2;
            }
        }

        previous.swap(current);
    }

    std::vector<std::string> aligned(profile_rows + 1);
    std::size_t i = rows - 1, j = columns - 1;

    while (i || j) {
        const auto move = direction[i * columns + j];

        if (move == 0 || move == 1) {
            --i;
            for (std::size_t row = 0; row < profile_rows; ++row) {
                aligned[row].push_back(profile[row][i]);
            }
        } else {
            for (std::size_t row = 0; row < profile_rows; ++row) {
                aligned[row].push_back('-');
            }
        }

        if (move == 0 || move == 2) aligned.back().push_back(sequence[--j]);
        else aligned.back().push_back('-');
    }

    for (auto& row : aligned) {
        std::reverse(row.begin(), row.end());
    }

    return aligned;
}

}

// Feng e Doolittle (1987), J. Mol. Evol. 25:351-360, PMID 3118049:
// alinhamento progressivo; neste projeto, a ordem é a dos acessos de montagem.
std::vector<std::string> align_progressively(const std::vector<std::string>& proteins) {
    if (proteins.empty())
        throw std::invalid_argument("No protein sequences to align");

    for (const auto& protein : proteins) {
        if (protein.empty()) throw std::invalid_argument("Empty protein sequence");
    }

    std::vector<std::string> aligned{proteins.front()};

    for (std::size_t index = 1; index < proteins.size(); ++index) {
        aligned = add_sequence(aligned, proteins[index]);
    }

    return aligned;
}

// Suyama, Torrents e Bork (2006), Nucleic Acids Res. 34:W609-W612,
// PMID 16845082: projeção do alinhamento proteico sobre os códons.
std::vector<std::string> project_multiple_codons(const std::vector<std::string>& aligned_proteins,
                                                 const std::vector<std::string>& coding_sequences) {
    if (aligned_proteins.empty() || aligned_proteins.size() != coding_sequences.size()) {
        throw std::invalid_argument("Protein and CDS counts differ");
    }

    const std::size_t columns = aligned_proteins.front().size();

    std::vector<std::string> projected;
    projected.reserve(aligned_proteins.size());

    for (std::size_t row = 0; row < aligned_proteins.size(); ++row) {
        if (aligned_proteins[row].size() != columns ||
            coding_sequences[row].size() != 3 * (residue_count(aligned_proteins[row]) + 1)) {
            throw std::invalid_argument("CDS length does not match protein alignment");
        }

        std::string nucleotides;
        nucleotides.reserve(3 * (columns + 1));

        std::size_t offset = 0;

        for (char amino_acid : aligned_proteins[row]) {
            if (amino_acid == '-') { 
                nucleotides += "---";
            } else {
                nucleotides += coding_sequences[row].substr(offset, 3);
                offset += 3;
            }
        }

        nucleotides += coding_sequences[row].substr(offset, 3);
        projected.push_back(std::move(nucleotides));
    }

    return projected;
}

ColumnSummary summarize_columns(const std::vector<std::string>& aligned_sequences) {
    if (aligned_sequences.empty()) throw std::invalid_argument("No aligned sequences");

    const std::size_t columns = aligned_sequences.front().size();
    ColumnSummary summary;

    for (const auto& row : aligned_sequences) {
        if (row.size() != columns)
            throw std::invalid_argument("Alignment rows have different lengths");
    }

    for (std::size_t column = 0; column < columns; ++column) {
        bool gap = false, same = true;

        for (const auto& row : aligned_sequences) {
            gap |= row[column] == '-';
            same &= row[column] == aligned_sequences.front()[column];
        }

        if (gap) ++summary.with_gap;
        else if (same) ++summary.conserved;
        else ++summary.substituted;
    }

    return summary;
}

// Zhang e Kumar (1997), Mol. Biol. Evol. 14:527-536, pp. 527-528:
// estas contagens são descritivas; convergência exige filogenia e ancestrais inferidos.
CodonEffects compare_codons(const Alignment& aligned_proteins,
                            const std::string& first_cds, const std::string& second_cds) {
    if (aligned_proteins.first.size() != aligned_proteins.second.size() ||
        first_cds.size() != 3 * (residue_count(aligned_proteins.first) + 1) ||
        second_cds.size() != 3 * (residue_count(aligned_proteins.second) + 1)) {
        throw std::invalid_argument("CDS length does not match paired proteins");
    }

    CodonEffects result;
    std::size_t first_offset = 0, second_offset = 0;

    for (std::size_t column = 0; column < aligned_proteins.first.size(); ++column) {
        const char amino_first = aligned_proteins.first[column];
        const char amino_second = aligned_proteins.second[column];

        if (amino_first == '-' || amino_second == '-') {
            ++result.gap_codons;

            if (amino_first != '-')
                first_offset += 3;

            if (amino_second != '-')
                second_offset += 3;

            continue;
        }

        const auto codon_first = first_cds.substr(first_offset, 3);
        const auto codon_second = second_cds.substr(second_offset, 3);

        for (std::size_t base = 0; base < 3; ++base) {
            result.different_bases += codon_first[base] != codon_second[base];
        }

        if (amino_first != amino_second) ++result.nonsynonymous_codons;
        else if (codon_first != codon_second) ++result.synonymous_codons;

        first_offset += 3;
        second_offset += 3;
    }

    result.different_stop = first_cds.substr(first_offset, 3) != second_cds.substr(second_offset, 3);

    return result;
}
