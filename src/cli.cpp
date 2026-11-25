#include "cli.hpp"

#include <ostream>
#include <stdexcept>

ProgramOptions parse_arguments(int argc, char* argv[], std::size_t available_threads) {
    ProgramOptions options;

    options.requested_threads = available_threads;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "--protein" && i + 1 < argc) {
            options.gene = argv[++i];

            if (options.gene != "recA" && options.gene != "trhO" && options.gene != "rpoB") {
                throw std::invalid_argument("--protein deve ser recA, trhO ou rpoB");
            }
        } else if (arg == "--input" && i + 1 < argc) {
            options.input_directory = argv[++i];
        } else if (arg == "--output-dir" && i + 1 < argc) {
            options.output_directory = argv[++i];
        } else if (arg == "--mode" && i + 1 < argc) {
            const std::string mode = argv[++i];

            if (mode == "protein") options.protein_mode = true;
            else if (mode == "cds") options.protein_mode = false;
            else throw std::invalid_argument("O modo deve ser protein ou cds");
        } else if (arg == "--threads" && i + 1 < argc) {
            const std::string value = argv[++i];
            std::size_t consumed = 0;

            if (value.empty() || value[0] == '-') {
                throw std::invalid_argument("Quantidade de threads inválida: " + value);
            }

            try {
                options.requested_threads = std::stoul(value, &consumed);
            } catch (const std::exception&) {
                throw std::invalid_argument("Quantidade de threads inválida: " + value);
            }

            if (consumed != value.size() || options.requested_threads == 0)
                throw std::invalid_argument("Quantidade de threads inválida: " + value);
        } else if (arg == "--benchmark") {
            options.benchmark = true;
        } else if (arg == "--help") {
            options.help = true;
            return options;
        } else {
            throw std::invalid_argument("Argumento desconhecido ou incompleto: " + arg);
        }
    }

    if (options.input_directory.empty())
        options.input_directory = std::filesystem::path("resources") / options.gene;

    if (options.output_directory.empty())
        options.output_directory = std::filesystem::path("results") / options.gene;

    return options;
}

void print_usage(std::ostream& output) {
    output << "Uso: gene-comparison [--protein recA|trhO|rpoB] "
              "[--mode protein|cds] [--threads N] [--benchmark] "
              "[--input DIRETÓRIO] [--output-dir DIRETÓRIO]\n";
}
