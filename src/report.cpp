#include "report.hpp"

#include "multiple_alignment.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace fs = std::filesystem;

void write_reports(const fs::path& directory, const std::string& gene,
                   const GeneSequences& data, const std::vector<PairResult>& pairs,
                   const std::vector<std::string>& aligned_proteins,
                   const std::vector<std::string>& aligned_cds) {
    fs::create_directories(directory);

    write_alignment(directory / "protein_alignment.faa", data.proteins, aligned_proteins, gene);
    write_alignment(directory / "cds_alignment.fna", data.coding_sequences, aligned_cds, gene);

    std::ofstream pair_file(directory / "pairs.tsv");

    if (!pair_file)
        throw std::runtime_error("Não foi possível criar pairs.tsv");

    pair_file << "montagem_1\tmontagem_2\tproteína_1\tproteína_2\tcoincidências_proteína\t"
                 "colunas_proteína\tidentidade_proteína_percentual\tlacunas_proteína\tcoincidências_cds\t"
                 "colunas_cds\tidentidade_cds_percentual\tlacunas_cds\n";

    std::ofstream codon_file(directory / "codon_effects.tsv");

    if (!codon_file)
        throw std::runtime_error("Não foi possível criar codon_effects.tsv");

    codon_file << "montagem_1\tmontagem_2\tbases_diferentes\tcódons_sinônimos\t"
                  "códons_não_sinônimos\tcódons_com_lacuna\tparada_diferente\n";

    for (const auto& pair : pairs) {
        const auto i = pair.first_index, j = pair.second_index;
        const auto& protein = pair.alignment;
        const auto cds = project_codons(protein, data.coding_sequences[i].sequence, data.coding_sequences[j].sequence);
        const auto effects = compare_codons(protein, data.coding_sequences[i].sequence, data.coding_sequences[j].sequence);

        pair_file << data.proteins[i].assembly << '\t' << data.proteins[j].assembly << '\t'
                  << data.proteins[i].accession << '\t' << data.proteins[j].accession << '\t'
                  << protein.matches << '\t' << protein.first.size() << '\t'
                  << std::fixed << std::setprecision(2) << protein.identity_percent() << '\t'
                  << protein.gaps << '\t' << cds.matches << '\t' << cds.first.size() << '\t'
                  << cds.identity_percent() << '\t' << cds.gaps << '\n';

        codon_file << data.proteins[i].assembly << '\t' << data.proteins[j].assembly << '\t'
                   << effects.different_bases << '\t' << effects.synonymous_codons << '\t'
                   << effects.nonsynonymous_codons << '\t' << effects.gap_codons << '\t'
                   << (effects.different_stop ? "sim" : "não") << '\n';
    }

    const auto protein_summary = summarize_columns(aligned_proteins);
    const auto cds_summary = summarize_columns(aligned_cds);

    std::ofstream summary_file(directory / "summary.tsv");

    if (!summary_file)
        throw std::runtime_error("Não foi possível criar summary.tsv");

    summary_file << "gene\tsequências\tcolunas_proteína\tproteína_conservada\tproteína_substituída\t"
                    "proteína_com_lacuna\tcolunas_cds\tcds_conservada\tcds_substituída\tcds_com_lacuna\n"
                 << gene << '\t' << data.proteins.size() << '\t' << aligned_proteins.front().size() << '\t'
                 << protein_summary.conserved << '\t' << protein_summary.substituted << '\t'
                 << protein_summary.with_gap << '\t' << aligned_cds.front().size() << '\t'
                 << cds_summary.conserved << '\t' << cds_summary.substituted << '\t'
                 << cds_summary.with_gap << '\n';
}

void print_console_report(const std::string& gene, bool protein_mode,
                          std::size_t threads, const fs::path& directory,
                          const GeneSequences& data, const std::vector<PairResult>& pairs,
                          const std::vector<std::string>& aligned_proteins) {
    std::cout << "gene=" << gene << " modo=" << (protein_mode ? "proteína" : "cds")
              << " sequências=" << data.proteins.size()
              << " pares_únicos=" << pairs.size() << " threads=" << threads << '\n';

    std::cout << "montagem_1\tmontagem_2\tacesso_1\tacesso_2\tcoincidências\tcolunas\t"
                 "identidade_percentual\tlacunas\tpontuação\tbases_diferentes\tcódons_sinônimos\t"
                 "códons_não_sinônimos\tcódons_com_lacuna\n";

    for (const auto& result : pairs) {
        const auto i = result.first_index, j = result.second_index;
        const auto& a = protein_mode ? data.proteins[i] : data.coding_sequences[i];
        const auto& b = protein_mode ? data.proteins[j] : data.coding_sequences[j];
        const auto alignment = protein_mode ? result.alignment : project_codons(result.alignment, data.cds_sequences[i], data.cds_sequences[j]);
        const auto effects = compare_codons(result.alignment, data.cds_sequences[i], data.cds_sequences[j]);

        std::cout << a.assembly << '\t' << b.assembly << '\t' << a.accession << '\t'
                  << b.accession << '\t' << alignment.matches << '\t' << alignment.first.size() << '\t'
                  << std::fixed << std::setprecision(2) << alignment.identity_percent() << '\t'
                  << alignment.gaps << '\t' << alignment.score << '\t' << effects.different_bases << '\t'
                  << effects.synonymous_codons << '\t' << effects.nonsynonymous_codons << '\t'
                  << effects.gap_codons << '\n';
    }

    const auto summary = summarize_columns(aligned_proteins);

    std::cout << "colunas_alinhamento_múltiplo=" << aligned_proteins.front().size()
              << " conservadas=" << summary.conserved << " substituídas=" << summary.substituted
              << " com_lacuna=" << summary.with_gap << '\n';

    std::cout << "diretório_resultados=" << directory.string() << '\n';
}
