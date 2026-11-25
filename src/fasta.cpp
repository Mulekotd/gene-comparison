#include "fasta.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>

namespace fs = std::filesystem;

namespace {

FastaRecord read_fasta(const fs::path& path, bool protein, const std::string& gene) {
    std::ifstream input(path);

    if (!input)
        throw std::runtime_error("Não foi possível abrir " + path.string());

    std::string header, line, sequence;

    if (!std::getline(input, header) || header.empty() || header[0] != '>') {
        throw std::runtime_error("Cabeçalho FASTA ausente em " + path.string());
    }

    while (std::getline(input, line)) {
        if (!line.empty() && line[0] == '>') {
            throw std::runtime_error("Era esperada uma sequência por arquivo: " + path.string());
        }

        for (unsigned char ch : line) {
            if (std::isspace(ch)) continue;

            const char letter = static_cast<char>(std::toupper(ch));
            const std::string alphabet = protein ? "ACDEFGHIKLMNPQRSTVWY" : "ACGT";

            if (alphabet.find(letter) == std::string::npos) {
                throw std::runtime_error("Caractere inválido na sequência em " + path.string());
            }

            sequence.push_back(letter);
        }
    }

    if (sequence.empty())
        throw std::runtime_error("Sequência FASTA vazia em " + path.string());

    const std::size_t end = header.find_first_of(" \t", 1);
    const std::string accession = header.substr(1, end == std::string::npos ? end : end - 1);

    std::string assembly = path.stem().string();
    const std::string suffix = "_" + gene;

    if (assembly.size() <= suffix.size() ||
        assembly.substr(assembly.size() - suffix.size()) != suffix) {
        throw std::runtime_error("Nome FASTA incompatível com " + gene + ": " + path.string());
    }

    assembly.resize(assembly.size() - suffix.size());
    return {assembly, accession, sequence};
}

std::vector<FastaRecord> load_records(const fs::path& directory, bool protein, const std::string& gene) {
    if (!fs::is_directory(directory))
        throw std::runtime_error("O caminho não é um diretório: " + directory.string());

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

}

GeneSequences load_gene_sequences(const fs::path& directory, const std::string& gene) {
    GeneSequences data;

    data.proteins = load_records(directory, true, gene);
    data.coding_sequences = load_records(directory, false, gene);

    for (std::size_t i = 0; i < data.proteins.size(); ++i) {
        if (data.coding_sequences[i].assembly != data.proteins[i].assembly) {
            throw std::runtime_error("Os FASTA de proteínas e CDS não correspondem");
        }

        data.protein_sequences.push_back(data.proteins[i].sequence);
        data.cds_sequences.push_back(data.coding_sequences[i].sequence);
    }

    return data;
}

void write_alignment(const fs::path& path, const std::vector<FastaRecord>& records,
                     const std::vector<std::string>& aligned, const std::string& gene) {
    std::ofstream output(path);

    if (!output)
        throw std::runtime_error("Não foi possível criar " + path.string());

    for (std::size_t row = 0; row < records.size(); ++row) {
        output << '>' << records[row].assembly << " acesso=" << records[row].accession
               << " gene=" << gene << '\n';

        for (std::size_t offset = 0; offset < aligned[row].size(); offset += 70) {
            output << aligned[row].substr(offset, 70) << '\n';
        }
    }
}
