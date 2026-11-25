#!/usr/bin/env python3
"""Extrai as CDS e proteínas do painel a partir dos genomas locais."""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

BASES = "TCAG"
COMPLEMENT = str.maketrans("ACGT", "TGCA")

AMINO_ACIDS = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)

CODON_TABLE = {
    a + b + c: amino
    for a, block in zip(BASES, range(0, 64, 16))
    for b, row in zip(BASES, range(block, block + 16, 4))
    for c, amino in zip(BASES, AMINO_ACIDS[row : row + 4])
}


def load_cds(row: dict[str, str]) -> tuple[str, str]:
    assembly = row["assembly"]
    files = list((ROOT / "resources" / "genomes").glob(f"{assembly}_*_genomic.fna"))

    if len(files) != 1:
        raise ValueError(f"Esperado um genoma para {assembly}; encontrados {len(files)}")

    with files[0].open(encoding="ascii") as handle:
        header = handle.readline()

        if not header.startswith(f'>{row["chromosome"]} '):
            raise ValueError(f"Cromossomo incorreto em {files[0]}")

        genome = "".join(line.strip().upper() for line in handle)

    start, end = int(row["start"]), int(row["end"])
    cds = genome[start - 1 : end]

    if row["strand"] == "-":
        cds = cds.translate(COMPLEMENT)[::-1]

    if not cds or len(cds) % 3 or set(cds) - set("ACGT"):
        raise ValueError(f"CDS inválida: {row['gene']} em {assembly}")

    protein = "".join(CODON_TABLE[cds[i : i + 3]] for i in range(0, len(cds), 3))

    if protein[0] != "M" or protein[-1] != "*" or "*" in protein[:-1]:
        raise ValueError(f"Tradução inválida: {row['gene']} em {assembly}")

    protein = protein[:-1]

    if len(protein) != int(row["amino_acids"]):
        raise ValueError(f"Comprimento inesperado: {row['gene']} em {assembly}")

    return cds, protein


def write_fasta(path: Path, header: str, sequence: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{header}\n")
        for offset in range(0, len(sequence), 70):
            handle.write(sequence[offset : offset + 70] + "\n")


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    with (ROOT / "resources" / "proteins.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    seen = set()

    for row in rows:
        gene, assembly = row["gene"], row["assembly"]

        if (gene, assembly) in seen:
            raise ValueError(f"Entrada duplicada: {gene} em {assembly}")

        seen.add((gene, assembly))
        cds, protein = load_cds(row)

        out = ROOT / "resources" / gene
        out.mkdir(parents=True, exist_ok=True)

        location = f'{row["chromosome"]}:{row["start"]}-{row["end"]}({row["strand"]})'
        metadata = f"montagem={assembly} gene={gene} id_proteína={row['protein_id']}"

        write_fasta(out / f"{assembly}_{gene}.fna", f"{location} {metadata}", cds)
        write_fasta(out / f"{assembly}_{gene}.faa", f"{row['protein_id']} {metadata} localização={location}", protein)

        print(f"{gene} {assembly}: {len(cds)} bases, {len(protein)} aminoácidos")

    if len(seen) != 12:
        raise ValueError(f"Esperadas 12 sequências no painel; encontradas {len(seen)}")


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
