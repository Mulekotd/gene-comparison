#!/usr/bin/env python3
"""Confere as saídas dos três genes com os FASTA e o cálculo independente."""

import csv
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROGRAM = ROOT / "build" / "gene-comparison"


def fasta(path: Path) -> dict[str, str]:
    records = {}
    name = None

    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        else:
            records[name] += line.strip()

    return records


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    check = subprocess.run(["python3", "scripts/analyze_proteins.py"], cwd=ROOT,
                           text=True, capture_output=True, check=True)
    independent = {}

    for row in csv.DictReader(check.stdout.splitlines()):
        key = (row["gene"], frozenset((row["referencia"], row["comparado"])))
        independent[key] = row

    with tempfile.TemporaryDirectory() as temporary:
        temporary = Path(temporary)

        for gene in ("recA", "trhO", "rpoB"):
            output = temporary / gene
            run = subprocess.run([str(PROGRAM), "--protein", gene, "--benchmark",
                                  "--output-dir", str(output)], cwd=ROOT,
                                 text=True, capture_output=True, check=True)

            assert "unique_pairs=6" in run.stdout
            assert "results_identical=yes" in run.stderr

            pairs = tsv(output / "pairs.tsv")
            effects = tsv(output / "codon_effects.tsv")
            summary = tsv(output / "summary.tsv")

            assert len(pairs) == len(effects) == 6
            assert len(summary) == 1
            assert len({frozenset((row["assembly_1"], row["assembly_2"])) for row in pairs}) == 6

            for row in effects:
                key = (gene, frozenset((row["assembly_1"], row["assembly_2"])))

                if key in independent:
                    expected = independent[key]

                    assert row["different_bases"] == expected["bases_diferentes"]
                    assert row["synonymous_codons"] == expected["codons_sinonimos"]
                    assert row["nonsynonymous_codons"] == expected["codons_nao_sinonimos"]

            protein_alignment = fasta(output / "protein_alignment.faa")
            cds_alignment = fasta(output / "cds_alignment.fna")

            assert len(protein_alignment) == len(cds_alignment) == 4

            protein_columns = {len(row) for row in protein_alignment.values()}
            cds_columns = {len(row) for row in cds_alignment.values()}

            assert len(protein_columns) == len(cds_columns) == 1
            assert next(iter(cds_columns)) == 3 * (next(iter(protein_columns)) + 1)

            for assembly, aligned in protein_alignment.items():
                original_protein = fasta(ROOT / "resources" / gene / f"{assembly}_{gene}.faa")
                original_cds = fasta(ROOT / "resources" / gene / f"{assembly}_{gene}.fna")

                assert aligned.replace("-", "") == next(iter(original_protein.values()))

                dna = cds_alignment[assembly]

                assert dna.replace("-", "") == next(iter(original_cds.values()))

                for column, amino_acid in enumerate(aligned):
                    codon = dna[3 * column : 3 * column + 3]
                    assert (codon == "---") == (amino_acid == "-")

            print(f"{gene}: seis pares e alinhamentos múltiplos validados")

        default = subprocess.run([str(PROGRAM), "--output-dir", str(temporary / "default")],
                                 cwd=ROOT, text=True, capture_output=True, check=True)

        assert "gene=recA" in default.stdout

        cds_mode = subprocess.run([str(PROGRAM), "--protein", "rpoB", "--mode", "cds",
                                   "--output-dir", str(temporary / "cds")], cwd=ROOT,
                                  text=True, capture_output=True, check=True)

        assert "gene=rpoB mode=cds" in cds_mode.stdout
        assert len(tsv(temporary / "cds" / "pairs.tsv")) == 6

        invalid = subprocess.run([str(PROGRAM), "--protein", "invalid"], cwd=ROOT,
                                 text=True, capture_output=True)

        assert invalid.returncode != 0

        invalid_threads = subprocess.run([str(PROGRAM), "--threads", "-1"], cwd=ROOT,
                                         text=True, capture_output=True)

        assert invalid_threads.returncode != 0


if __name__ == "__main__":
    main()
