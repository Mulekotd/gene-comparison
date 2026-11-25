# Resultados iniciais do painel

O programa comparou as CDS e proteínas de `recA`, `trhO` e `rpoB` nas quatro linhagens do projeto. Em cada gene, calculou os seis pares únicos com threads e construiu um alinhamento múltiplo progressivo das proteínas, depois projetado sobre as CDS. A ordem de inclusão das linhagens no alinhamento múltiplo é a ordem alfabética dos acessos de montagem; esse método é uma heurística e pode depender da ordem das sequências.

| Gene | Aminoácidos por proteína | Identidade entre pares | Colunas proteicas conservadas | Colunas com substituição | Colunas com lacuna | Bases conservadas na CDS alinhada |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `recA` | 335 | 94,93–97,91% | 312 | 23 | 0 | 693/1.008 |
| `trhO` | 343 | 86,30–93,59% | 282 | 61 | 0 | 651/1.032 |
| `rpoB` | 1.269 | 96,45–98,11% | 1.202 | 67 | 0 | 2.815/3.810 |

Neste painel, `trhO` apresenta a menor identidade proteica entre os pares e a maior proporção de colunas proteicas variáveis. Nenhuma das três proteínas apresentou lacunas nos alinhamentos das quatro linhagens. O DNA varia mais que a proteína porque códons diferentes podem codificar o mesmo aminoácido. Por exemplo, na comparação de `recA` entre _M. aestuarii_ e a linhagem 1_2014MBL_MicDiv, há 213 bases diferentes nas regiões codificadoras: 173 códons diferentes preservam o aminoácido e 13 trocam o aminoácido. As unidades dessas contagens diferem; um códon pode conter mais de uma base diferente.

Os valores de identidade são proporções de posições iguais nos alinhamentos entre pares, com lacunas incluídas no denominador. As colunas conservadas e variáveis vêm do alinhamento múltiplo. O códon de parada aparece na CDS e não na proteína; diferenças entre códons de parada são relatadas separadamente em `codon_effects.tsv`.

Estas são diferenças observadas entre sequências atuais. Não mostram em qual linhagem uma alteração ocorreu, nem demonstram mudança de função da proteína. A [tabela de proteínas](proteinas.md) apresenta os acessos RefSeq e exemplos de códons; a [avaliação de ortologia](ortologia.md) descreve as evidências para os três genes.

Arquivos produzidos para cada gene em `results/<gene>/`:

- `protein_alignment.faa` e `cds_alignment.fna`: alinhamentos múltiplos;
- `pairs.tsv`: identidades e lacunas nos seis pares;
- `codon_effects.tsv`: bases diferentes, códons sinônimos, não sinônimos e com lacuna;
- `summary.tsv`: colunas conservadas, substituídas e com lacuna.
