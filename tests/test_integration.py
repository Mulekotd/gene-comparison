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

    for line in path.read_text(encoding="utf-8").splitlines():
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
        key = (row["gene"], frozenset((row["referência"], row["comparado"])))
        independent[key] = row

    with tempfile.TemporaryDirectory() as temporary:
        temporary = Path(temporary)

        for gene in ("recA", "trhO", "rpoB"):
            output = temporary / gene
            run = subprocess.run([str(PROGRAM), "--protein", gene, "--benchmark",
                                  "--output-dir", str(output)], cwd=ROOT,
                                 text=True, capture_output=True, check=True)

            assert "pares_únicos=6" in run.stdout
            assert "resultados_idênticos=sim" in run.stderr

            pairs = tsv(output / "pairs.tsv")
            effects = tsv(output / "codon_effects.tsv")
            summary = tsv(output / "summary.tsv")

            assert len(pairs) == len(effects) == 6
            assert len(summary) == 1
            assert "proteína_1" in pairs[0]
            assert "códons_sinônimos" in effects[0]
            assert "sequências" in summary[0]
            assert {row["parada_diferente"] for row in effects} <= {"sim", "não"}
            assert len({frozenset((row["montagem_1"], row["montagem_2"])) for row in pairs}) == 6

            for row in effects:
                key = (gene, frozenset((row["montagem_1"], row["montagem_2"])))

                if key in independent:
                    expected = independent[key]

                    assert row["bases_diferentes"] == expected["bases_diferentes"]
                    assert row["códons_sinônimos"] == expected["códons_sinônimos"]
                    assert row["códons_não_sinônimos"] == expected["códons_não_sinônimos"]

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

        assert "gene=rpoB modo=cds" in cds_mode.stdout
        assert len(tsv(temporary / "cds" / "pairs.tsv")) == 6

        invalid = subprocess.run([str(PROGRAM), "--protein", "invalid"], cwd=ROOT,
                                 text=True, capture_output=True)

        assert invalid.returncode != 0
        assert "erro: --protein deve ser recA, trhO ou rpoB" in invalid.stderr

        invalid_threads = subprocess.run([str(PROGRAM), "--threads", "-1"], cwd=ROOT,
                                         text=True, capture_output=True)

        assert invalid_threads.returncode != 0
        assert "erro: Quantidade de threads inválida: -1" in invalid_threads.stderr

        overflow_threads = subprocess.run([str(PROGRAM), "--threads", "999999999999999999999999"],
                                          cwd=ROOT, text=True, capture_output=True)

        assert overflow_threads.returncode != 0
        assert "erro: Quantidade de threads inválida" in overflow_threads.stderr

        invalid_mode = subprocess.run([str(PROGRAM), "--mode", "desconhecido"],
                                      cwd=ROOT, text=True, capture_output=True)

        assert invalid_mode.returncode != 0
        assert "erro: O modo deve ser protein ou cds" in invalid_mode.stderr


if __name__ == "__main__":
    main()
