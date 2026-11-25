#pragma once

#include "alignment.hpp"
#include "fasta.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

void write_reports(const std::filesystem::path& directory, const std::string& gene,
                   const GeneSequences& data, const std::vector<PairResult>& pairs,
                   const std::vector<std::string>& aligned_proteins,
                   const std::vector<std::string>& aligned_cds);

void print_console_report(const std::string& gene, bool protein_mode,
                          std::size_t threads, const std::filesystem::path& directory,
                          const GeneSequences& data, const std::vector<PairResult>& pairs,
                          const std::vector<std::string>& aligned_proteins);
