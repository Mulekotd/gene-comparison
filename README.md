# Comparação de um gene em quatro espécies de _Maribacter_

Projeto de bioinformática do tema 4: comparar um gene ortólogo em quatro espécies do mesmo gênero e avaliar diferenças nas sequências de DNA e nas proteínas correspondentes. A entrega será o código em C++ e um relatório em formato de artigo, conforme [o modelo SBC](docs/template-latex/sbc-template.tex).

## Espécies escolhidas

Selecionamos uma linhagem por espécie a partir de [docs/samples.csv](docs/samples.csv). Todas constam como espécies validadas e isolados com genoma completo em um contig.

| Espécie             | Acesso do genoma  | Linhagem         | Completude | Contaminação |
| ------------------- | ----------------- | ---------------- | ---------: | -----------: |
| _M. aestuarii_      | `GCF_027474845.2` | JCM18631         |     97,57% |        0,07% |
| _M. cobaltidurans_  | `GCF_002269385.1` | B1               |     99,89% |        0,32% |
| _M. dokdonensis_    | `GCF_001865655.1` | 1_2014MBL_MicDiv |     99,99% |        0,02% |
| _M. hydrothermalis_ | `GCF_001913155.1` | T28              |     99,99% |        0,27% |

A planilha contém apenas metadados. A [avaliação de ortologia de `recA`](docs/ortologia-reca.md) foi registrada separadamente; as sequências de entrada ainda precisam ser reunidas em FASTA.

## Gene e dados

O **gene escolhido é `recA`**, envolvido na recombinação e no reparo do DNA. Sua CDS tem tamanho adequado para iniciar o alinhamento por programação dinâmica. A avaliação das quatro linhagens encontrou uma CDS íntegra e uma cópia correspondente em cada montagem.

Para cada linhagem, extrair a CDS e a proteína de `recA` em **FASTA** a partir do registro do genoma no [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/). FASTA será o formato de entrada do programa. Os acessos, as posições, as orientações e a verificação das traduções estão na [avaliação de ortologia](docs/ortologia-reca.md).

As próximas tarefas estão em [docs/desenvolvimento.md](docs/desenvolvimento.md).
