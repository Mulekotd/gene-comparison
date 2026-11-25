#include "alignment.hpp"
#include "cli.hpp"
#include "fasta.hpp"
#include "multiple_alignment.hpp"
#include "report.hpp"
#include "utils.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <thread>

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    try {
        const std::size_t available_threads = std::max(1u, std::thread::hardware_concurrency());
        const auto options = parse_arguments(argc, argv, available_threads);

        if (options.help) {
            print_usage(std::cout);
            return 0;
        }

        const auto data = load_gene_sequences(options.input_directory, options.gene);
        const std::size_t pair_count = unique_pairs(data.protein_sequences.size()).size();
        const std::size_t threads = std::min({options.requested_threads, available_threads, pair_count});

        const auto start = std::chrono::steady_clock::now();
        const auto results = align_pairs(data.protein_sequences, threads);
        const auto parallel_end = std::chrono::steady_clock::now();
        const auto multiple_proteins = align_progressively(data.protein_sequences);
        const auto multiple_cds = project_multiple_codons(multiple_proteins, data.cds_sequences);

        write_reports(options.output_directory, options.gene, data, results, multiple_proteins, multiple_cds);
        print_console_report(options.gene, options.protein_mode, threads,
                             options.output_directory, data, results, multiple_proteins);

        if (options.benchmark) {
            const auto sequential_start = std::chrono::steady_clock::now();
            const auto serial = align_pairs(data.protein_sequences, 1);
            const auto sequential_end = std::chrono::steady_clock::now();

            if (!same_results(results, serial))
                throw std::runtime_error("Os resultados paralelos e sequenciais são diferentes");

            const auto parallel_ms = std::chrono::duration<double, std::milli>(parallel_end - start).count();
            const auto serial_ms = std::chrono::duration<double, std::milli>(sequential_end - sequential_start).count();

            std::cerr << std::fixed << std::setprecision(3)
                      << "tempo_pares_paralelo_ms=" << parallel_ms << " tempo_pares_sequencial_ms=" << serial_ms
                      << " resultados_idênticos=sim\n";
        }

        return 0;
    } catch (const fs::filesystem_error& error) {
        std::cerr << "erro: Falha no sistema de arquivos";
        if (!error.path1().empty()) std::cerr << " em " << error.path1().string();
        if (!error.path2().empty()) std::cerr << " e " << error.path2().string();
        std::cerr << " (código " << error.code().value() << ")\n";
        return 1;
    } catch (const std::exception& error) {
        std::cerr << "erro: " << error.what() << '\n';
        return 1;
    }
}
