#!/usr/bin/env python3
"""Conta diferenças de códons nas proteínas do painel, usando os genomas locais."""

import csv
import sys
from collections import defaultdict
from pathlib import Path

from extract_proteins import load_cds


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = "GCF_027474845.2"


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    by_gene = defaultdict(dict)

    with (ROOT / "resources" / "proteins.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            gene, assembly = row["gene"], row["assembly"]

            if assembly in by_gene[gene]:
                raise ValueError(f"Entrada duplicada: {gene} em {assembly}")

            by_gene[gene][assembly] = load_cds(row)

    writer = csv.writer(sys.stdout)
    writer.writerow(("gene", "referência", "comparado", "bases_diferentes", "códons_sinônimos",
                     "códons_não_sinônimos", "aminoácidos_diferentes", "códon_parada_referência",
                     "códon_parada_comparado"))

    for gene, sequences in sorted(by_gene.items()):
        if len(sequences) != 4 or REFERENCE not in sequences:
            raise ValueError(f"Painel incompleto para {gene}")

        ref_cds, ref_protein = sequences[REFERENCE]

        for assembly, (cds, protein) in sorted(sequences.items()):
            if assembly == REFERENCE:
                continue

            if len(cds) != len(ref_cds) or len(protein) != len(ref_protein):
                raise ValueError(f"Alinhamento de códons necessário para {gene}: {assembly}")

            base_differences = sum(a != b for a, b in zip(ref_cds[:-3], cds[:-3]))

            synonymous = sum(
                ref_cds[3 * i : 3 * i + 3] != cds[3 * i : 3 * i + 3]
                and ref_protein[i] == protein[i]
                for i in range(len(protein))
            )

            amino_acid_differences = sum(a != b for a, b in zip(ref_protein, protein))

            writer.writerow((gene, REFERENCE, assembly, base_differences, synonymous,
                             amino_acid_differences, amino_acid_differences,
                             ref_cds[-3:], cds[-3:]))


if __name__ == "__main__":
    sys.stderr.reconfigure(encoding="utf-8")

    try:
        main()
    except ValueError as erro:
        print(f"erro: {erro}", file=sys.stderr)
        raise SystemExit(1) from None
    except KeyError as erro:
        print(f"erro: Campo ausente nos dados: {erro.args[0]}", file=sys.stderr)
        raise SystemExit(1) from None
    except OSError as erro:
        print(f"erro: Falha ao acessar o sistema de arquivos em {erro.filename} (código {erro.errno})", file=sys.stderr)
        raise SystemExit(1) from None
