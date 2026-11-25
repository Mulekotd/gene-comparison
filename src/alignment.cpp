#include "alignment.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>

// Yang (2006), Computational Molecular Evolution, cap. 1, pp. 3-4:
// a proporção observada de diferenças não é a distância evolutiva corrigida.
double Alignment::identity_percent() const {
    return first.empty() ? 0.0 : 100.0 * static_cast<double>(matches) / first.size();
}

// Needleman e Wunsch (1970), J. Mol. Biol. 48:443-453, PMID 5420325:
// alinhamento global de duas sequências por programação dinâmica.
Alignment align_global(const std::string& first, const std::string& second, int match, int mismatch, int gap) {
    if (first.empty() || second.empty()) {
        throw std::invalid_argument("Cannot align an empty sequence");
    }

    const std::size_t rows = first.size() + 1;
    const std::size_t columns = second.size() + 1;

    // Um byte por célula guarda o caminho; os scores precisam de apenas duas linhas.
    std::vector<std::uint8_t> direction(rows * columns, 0);
    std::vector<int> previous(columns), current(columns);

    for (std::size_t j = 1; j < columns; ++j) {
        previous[j] = previous[j - 1] + gap;
        direction[j] = 2;  // esquerda
    }

    for (std::size_t i = 1; i < rows; ++i) {
        current[0] = previous[0] + gap;
        direction[i * columns] = 1;  // acima

        for (std::size_t j = 1; j < columns; ++j) {
            const int diagonal = previous[j - 1] + (first[i - 1] == second[j - 1] ? match : mismatch);
            const int up = previous[j] + gap;
            const int left = current[j - 1] + gap;

            // Em caso de empate, prioriza a diagonal para reproduzir o alinhamento.
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

    Alignment result;
    result.score = previous.back();

    std::size_t i = first.size(), j = second.size();

    result.first.reserve(i + j);
    result.second.reserve(i + j);

    while (i || j) {
        const auto move = direction[i * columns + j];

        if (move == 0) {
            result.first.push_back(first[--i]);
            result.second.push_back(second[--j]);
            result.matches += result.first.back() == result.second.back();
        } else if (move == 1) {
            result.first.push_back(first[--i]);
            result.second.push_back('-');
            ++result.gaps;
        } else {
            result.first.push_back('-');
            result.second.push_back(second[--j]);
            ++result.gaps;
        }
    }

    std::reverse(result.first.begin(), result.first.end());
    std::reverse(result.second.begin(), result.second.end());

    return result;
}

// Suyama, Torrents e Bork (2006), Nucleic Acids Res. 34:W609-W612,
// PMID 16845082: correspondência entre resíduos alinhados e códons da CDS.
Alignment project_codons(const Alignment& protein_alignment,
                         const std::string& first_cds, const std::string& second_cds) {
    const auto count_residues = [](const std::string& aligned) {
        return static_cast<std::size_t>(std::count_if(aligned.begin(), aligned.end(),
                                                       [](char ch) { return ch != '-'; }));
    };

    if (protein_alignment.first.size() != protein_alignment.second.size() ||
        first_cds.size() != 3 * (count_residues(protein_alignment.first) + 1) ||
        second_cds.size() != 3 * (count_residues(protein_alignment.second) + 1)) {
        throw std::invalid_argument("CDS length does not match protein alignment plus stop codon");
    }

    Alignment result;
    std::size_t first_offset = 0, second_offset = 0;

    for (std::size_t column = 0; column < protein_alignment.first.size(); ++column) {
        if (protein_alignment.first[column] == '-') {
            result.first += "---";
        } else {
            result.first += first_cds.substr(first_offset, 3);
            first_offset += 3;
        }

        if (protein_alignment.second[column] == '-') {
            result.second += "---";
        } else {
            result.second += second_cds.substr(second_offset, 3);
            second_offset += 3;
        }
    }

    // Os códons de parada fazem parte das CDS, mas não dos FASTA de proteínas.
    result.first += first_cds.substr(first_offset, 3);
    result.second += second_cds.substr(second_offset, 3);

    for (std::size_t i = 0; i < result.first.size(); ++i) {
        const char a = result.first[i], b = result.second[i];

        if (a == '-' || b == '-') {
            ++result.gaps;
            result.score -= 3;
        } else if (a == b) {
            ++result.matches;
            result.score += 2;
        } else result.score -= 1;
    }

    return result;
}

std::vector<std::pair<std::size_t, std::size_t>> unique_pairs(std::size_t count) {
    std::vector<std::pair<std::size_t, std::size_t>> pairs;

    if (count < 2)
        return pairs;

    pairs.reserve(count * (count - 1) / 2);

    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count; ++j) {
            pairs.emplace_back(i, j);
        }
    }

    return pairs;
}

std::vector<PairResult> align_pairs(const std::vector<std::string>& sequences,
                                    std::size_t thread_count) {
    if (thread_count == 0) {
        throw std::invalid_argument("Thread count must be positive");
    }

    const auto pairs = unique_pairs(sequences.size());

    std::vector<PairResult> results(pairs.size());

    if (pairs.empty())
        return results;

    std::atomic<std::size_t> next_task{0};
    std::exception_ptr error;
    std::mutex error_mutex;

    const auto worker = [&] {
        while (true) {
            const std::size_t task = next_task.fetch_add(1);

            if (task >= pairs.size()) break;

            const auto [i, j] = pairs[task];

            try {
                results[task] = {i, j, align_global(sequences[i], sequences[j])};
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);

                if (!error) {
                    error = std::current_exception();
                }
            }
        }
    };

    const std::size_t workers = std::min(thread_count, pairs.size());
    std::vector<std::thread> threads;

    threads.reserve(workers > 1 ? workers - 1 : 0);

    for (std::size_t i = 1; i < workers; ++i) threads.emplace_back(worker);

    worker();

    for (auto& thread : threads) thread.join();

    if (error)
        std::rethrow_exception(error);

    return results;
}
