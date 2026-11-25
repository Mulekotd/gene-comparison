#include "alignment.hpp"
#include "multiple_alignment.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace fs = std::filesystem;

struct FastaRecord {
    std::string assembly;
    std::string accession;
    std::string sequence;
};

FastaRecord read_fasta(const fs::path& path, bool protein, const std::string& gene) {
    std::ifstream input(path);

    if (!input)
        throw std::runtime_error("Cannot open " + path.string());

    std::string header, line, sequence;

    if (!std::getline(input, header) || header.empty() || header[0] != '>') {
        throw std::runtime_error("Missing FASTA header in " + path.string());
    }

    while (std::getline(input, line)) {
        if (!line.empty() && line[0] == '>') {
            throw std::runtime_error("Expected one sequence per file: " + path.string());
        }

        for (unsigned char ch : line) {
            if (std::isspace(ch)) continue;

            const char letter = static_cast<char>(std::toupper(ch));
            const std::string alphabet = protein ? "ACDEFGHIKLMNPQRSTVWY" : "ACGT";

            if (alphabet.find(letter) == std::string::npos) {
                throw std::runtime_error("Invalid sequence character in " + path.string());
            }

            sequence.push_back(letter);
        }
    }

    if (sequence.empty())
        throw std::runtime_error("Empty FASTA sequence in " + path.string());

    const std::size_t end = header.find_first_of(" \t", 1);
    const std::string accession = header.substr(1, end == std::string::npos ? end : end - 1);

    std::string assembly = path.stem().string();

    const std::string suffix = "_" + gene;
    if (assembly.size() <= suffix.size() ||
        assembly.substr(assembly.size() - suffix.size()) != suffix) {
        throw std::runtime_error("Nome FASTA incompatível com " + gene + ": " + path.string());
    }

    assembly.resize(assembly.size() - suffix.size());

    return { assembly, accession, sequence };
}

std::vector<FastaRecord> load_records(const fs::path& directory, bool protein,
                                      const std::string& gene) {
    if (!fs::is_directory(directory))
        throw std::runtime_error("Not a directory: " + directory.string());

    std::vector<fs::path> files;

    const std::string extension = protein ? ".faa" : ".fna";

    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().extension() == extension) {
            files.push_back(entry.path());
        }
    }

    std::sort(files.begin(), files.end());

    if (files.size() != 4)
        throw std::runtime_error("Esperados quatro FASTA em " + directory.string());

    std::vector<FastaRecord> records;

    for (const auto& file : files) {
        records.push_back(read_fasta(file, protein, gene));
    }

    return records;
}

bool same_results(const std::vector<PairResult>& a, const std::vector<PairResult>& b) {
    if (a.size() != b.size())
        return false;

    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].first_index != b[i].first_index || a[i].second_index != b[i].second_index ||
            a[i].alignment.score != b[i].alignment.score ||
            a[i].alignment.first != b[i].alignment.first ||
            a[i].alignment.second != b[i].alignment.second) return false;
    }

    return true;
}

void write_alignment(const fs::path& path, const std::vector<FastaRecord>& records,
                     const std::vector<std::string>& aligned, const std::string& gene) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Não foi possível criar " + path.string());
    for (std::size_t row = 0; row < records.size(); ++row) {
        output << '>' << records[row].assembly << " accession=" << records[row].accession
               << " gene=" << gene << '\n';
        for (std::size_t offset = 0; offset < aligned[row].size(); offset += 70) {
            output << aligned[row].substr(offset, 70) << '\n';
        }
    }
}

void write_reports(const fs::path& directory, const std::string& gene,
                   const std::vector<FastaRecord>& proteins,
                   const std::vector<FastaRecord>& coding_sequences,
                   const std::vector<PairResult>& pairs,
                   const std::vector<std::string>& aligned_proteins,
                   const std::vector<std::string>& aligned_cds) {
    fs::create_directories(directory);
    write_alignment(directory / "protein_alignment.faa", proteins, aligned_proteins, gene);
    write_alignment(directory / "cds_alignment.fna", coding_sequences, aligned_cds, gene);

    std::ofstream pair_file(directory / "pairs.tsv");
    if (!pair_file) throw std::runtime_error("Não foi possível criar pairs.tsv");
    pair_file << "assembly_1\tassembly_2\tprotein_1\tprotein_2\tprotein_matches\t"
                 "protein_columns\tprotein_identity_percent\tprotein_gaps\tcds_matches\t"
                 "cds_columns\tcds_identity_percent\tcds_gaps\n";
    std::ofstream codon_file(directory / "codon_effects.tsv");
    if (!codon_file) throw std::runtime_error("Não foi possível criar codon_effects.tsv");
    codon_file << "assembly_1\tassembly_2\tdifferent_bases\tsynonymous_codons\t"
                  "nonsynonymous_codons\tgap_codons\tdifferent_stop\n";

    for (const auto& pair : pairs) {
        const auto i = pair.first_index, j = pair.second_index;
        const auto& protein = pair.alignment;
        const auto cds = project_codons(protein, coding_sequences[i].sequence,
                                       coding_sequences[j].sequence);
        const auto effects = compare_codons(protein, coding_sequences[i].sequence,
                                            coding_sequences[j].sequence);
        pair_file << proteins[i].assembly << '\t' << proteins[j].assembly << '\t'
                  << proteins[i].accession << '\t' << proteins[j].accession << '\t'
                  << protein.matches << '\t' << protein.first.size() << '\t'
                  << std::fixed << std::setprecision(2) << protein.identity_percent() << '\t'
                  << protein.gaps << '\t' << cds.matches << '\t' << cds.first.size() << '\t'
                  << cds.identity_percent() << '\t' << cds.gaps << '\n';
        codon_file << proteins[i].assembly << '\t' << proteins[j].assembly << '\t'
                   << effects.different_bases << '\t' << effects.synonymous_codons << '\t'
                   << effects.nonsynonymous_codons << '\t' << effects.gap_codons << '\t'
                   << (effects.different_stop ? "yes" : "no") << '\n';
    }

    const auto protein_summary = summarize_columns(aligned_proteins);
    const auto cds_summary = summarize_columns(aligned_cds);
    std::ofstream summary_file(directory / "summary.tsv");
    if (!summary_file) throw std::runtime_error("Não foi possível criar summary.tsv");
    summary_file << "gene\tsequences\tprotein_columns\tprotein_conserved\tprotein_substituted\t"
                    "protein_with_gap\tcds_columns\tcds_conserved\tcds_substituted\tcds_with_gap\n"
                 << gene << '\t' << proteins.size() << '\t' << aligned_proteins.front().size() << '\t'
                 << protein_summary.conserved << '\t' << protein_summary.substituted << '\t'
                 << protein_summary.with_gap << '\t' << aligned_cds.front().size() << '\t'
                 << cds_summary.conserved << '\t' << cds_summary.substituted << '\t'
                 << cds_summary.with_gap << '\n';
}

int main(int argc, char* argv[]) {
    try {
        std::string gene = "recA";
        fs::path input_directory, output_directory;
        bool protein = true, benchmark = false;

        const std::size_t available_threads = std::max(1u, std::thread::hardware_concurrency());
        std::size_t requested_threads = available_threads;

        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];

            if (arg == "--protein" && i + 1 < argc) {
                gene = argv[++i];
                if (gene != "recA" && gene != "trhO" && gene != "rpoB") {
                    throw std::invalid_argument("--protein deve ser recA, trhO ou rpoB");
                }
            } else if (arg == "--input" && i + 1 < argc) {
                input_directory = argv[++i];
            } else if (arg == "--output-dir" && i + 1 < argc) {
                output_directory = argv[++i];
            } else if (arg == "--mode" && i + 1 < argc) {
                const std::string mode = argv[++i];
                if (mode == "protein") protein = true;
                else if (mode == "cds") protein = false;
                else throw std::invalid_argument("Mode must be protein or cds");
            } else if (arg == "--threads" && i + 1 < argc) {
                const std::string value = argv[++i];
                std::size_t consumed = 0;
                if (value.empty() || value[0] == '-') {
                    throw std::invalid_argument("Invalid thread count: " + value);
                }
                requested_threads = std::stoul(value, &consumed);

                if (consumed != value.size() || requested_threads == 0)
                    throw std::invalid_argument("Invalid thread count: " + value);
            } else if (arg == "--benchmark") {
                benchmark = true;
            } else if (arg == "--help") {
                std::cout << "Uso: gene-comparison [--protein recA|trhO|rpoB] "
                             "[--mode protein|cds] [--threads N] [--benchmark] "
                             "[--input DIR] [--output-dir DIR]\n";
                return 0;
            } else {
                throw std::invalid_argument("Unknown or incomplete argument: " + arg);
            }
        }

        if (input_directory.empty()) input_directory = fs::path("resources") / gene;
        if (output_directory.empty()) output_directory = fs::path("results") / gene;
        const auto proteins = load_records(input_directory, true, gene);
        const auto coding_sequences = load_records(input_directory, false, gene);

        std::vector<std::string> sequences;
        std::vector<std::string> cds_sequences;

        for (std::size_t i = 0; i < proteins.size(); ++i) {
            if (coding_sequences[i].assembly != proteins[i].assembly) {
                throw std::runtime_error("Os FASTA de proteínas e CDS não correspondem");
            }
            sequences.push_back(proteins[i].sequence);
            cds_sequences.push_back(coding_sequences[i].sequence);
        }

        const std::size_t pair_count = unique_pairs(sequences.size()).size();
        const std::size_t threads = std::min({requested_threads, available_threads, pair_count});

        const auto start = std::chrono::steady_clock::now();
        const auto results = align_pairs(sequences, threads);
        const auto parallel_end = std::chrono::steady_clock::now();
        const auto multiple_proteins = align_progressively(sequences);
        const auto multiple_cds = project_multiple_codons(multiple_proteins, cds_sequences);
        write_reports(output_directory, gene, proteins, coding_sequences, results,
                      multiple_proteins, multiple_cds);

        std::cout << "gene=" << gene << " mode=" << (protein ? "protein" : "cds")
                  << " sequences=" << proteins.size()
                  << " unique_pairs=" << pair_count << " threads=" << threads << '\n';

        std::cout << "assembly_1\tassembly_2\taccession_1\taccession_2\tmatches\tcolumns\t"
                     "identity_percent\tgaps\tscore\tdifferent_bases\tsynonymous_codons\t"
                     "nonsynonymous_codons\tgap_codons\n";

        for (const auto& result : results) {
            const auto i = result.first_index, j = result.second_index;
            const auto& a = protein ? proteins[i] : coding_sequences[i];
            const auto& b = protein ? proteins[j] : coding_sequences[j];
            const auto aln = protein ? result.alignment
                                     : project_codons(result.alignment, cds_sequences[i], cds_sequences[j]);
            const auto effects = compare_codons(result.alignment, cds_sequences[i], cds_sequences[j]);

            std::cout << a.assembly << '\t' << b.assembly << '\t' << a.accession << '\t'
                      << b.accession << '\t' << aln.matches << '\t' << aln.first.size() << '\t'
                      << std::fixed << std::setprecision(2) << aln.identity_percent() << '\t'
                      << aln.gaps << '\t' << aln.score << '\t' << effects.different_bases << '\t'
                      << effects.synonymous_codons << '\t' << effects.nonsynonymous_codons << '\t'
                      << effects.gap_codons << '\n';
        }

        const auto summary = summarize_columns(multiple_proteins);
        std::cout << "multiple_alignment_columns=" << multiple_proteins.front().size()
                  << " conserved=" << summary.conserved << " substituted=" << summary.substituted
                  << " with_gap=" << summary.with_gap << '\n';
        std::cout << "results_dir=" << output_directory.string() << '\n';

        if (benchmark) {
            const auto sequential_start = std::chrono::steady_clock::now();
            const auto serial = align_pairs(sequences, 1);
            const auto sequential_end = std::chrono::steady_clock::now();

            if (!same_results(results, serial))
                throw std::runtime_error("Parallel and serial results differ");

            const auto parallel_ms = std::chrono::duration<double, std::milli>(parallel_end - start).count();
            const auto serial_ms = std::chrono::duration<double, std::milli>(sequential_end - sequential_start).count();

            std::cerr << std::fixed << std::setprecision(3)
                      << "pairwise_parallel_ms=" << parallel_ms << " pairwise_sequential_ms=" << serial_ms
                      << " results_identical=yes\n";
        }

        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
