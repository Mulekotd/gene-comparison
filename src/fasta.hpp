#pragma once

#include <filesystem>
#include <string>
#include <vector>

struct FastaRecord {
    std::string assembly;
    std::string accession;
    std::string sequence;
};

struct GeneSequences {
    std::vector<FastaRecord> proteins;
    std::vector<FastaRecord> coding_sequences;
    std::vector<std::string> protein_sequences;
    std::vector<std::string> cds_sequences;
};

GeneSequences load_gene_sequences(const std::filesystem::path& directory, const std::string& gene);
void write_alignment(const std::filesystem::path& path,
                     const std::vector<FastaRecord>& records,
                     const std::vector<std::string>& aligned,
                     const std::string& gene);
