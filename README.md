# Comparação de um gene em quatro espécies de _Maribacter_

Projeto de bioinformática do tema 4: comparar um gene ortólogo em quatro espécies do mesmo gênero e avaliar diferenças nas sequências de DNA e nas proteínas correspondentes. A entrega será o código em C++ e um relatório em formato de artigo, conforme o modelo SBC indicado pela professora.

## Espécies escolhidas

Selecionamos uma linhagem por espécie a partir de [docs/samples.csv](docs/samples.csv). Todas constam como espécies validadas e isolados com genoma completo em um contig.

| Espécie             | Acesso do genoma  | Linhagem         | Completude | Contaminação |
| ------------------- | ----------------- | ---------------- | ---------: | -----------: |
| _M. aestuarii_      | `GCF_027474845.2` | JCM18631         |     97,57% |        0,07% |
| _M. cobaltidurans_  | `GCF_002269385.1` | B1               |     99,89% |        0,32% |
| _M. dokdonensis_    | `GCF_001865655.1` | 1_2014MBL_MicDiv |     99,99% |        0,02% |
| _M. hydrothermalis_ | `GCF_001913155.1` | T28              |     99,99% |        0,27% |

A planilha contém apenas metadados. A [avaliação de ortologia dos três genes](docs/ortologia.md) foi registrada separadamente; os FASTA das CDS e proteínas analisadas estão em `resources/<gene>/`.

## Gene e dados

O **gene escolhido é `recA`**, envolvido na recombinação e no reparo do DNA. Sua CDS tem tamanho adequado para iniciar o alinhamento por programação dinâmica. A avaliação das quatro linhagens encontrou uma CDS íntegra e uma cópia correspondente em cada montagem.

As CDS e proteínas foram extraídas em **FASTA** dos genomas obtidos no [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/). FASTA é o formato de entrada do programa. Os acessos, as posições, as orientações e a verificação das traduções estão na [avaliação de ortologia](docs/ortologia.md) e em [proteins.csv](resources/proteins.csv).

As próximas tarefas estão em [docs/desenvolvimento.md](docs/desenvolvimento.md).

As passagens dos PDFs usados para orientar a interpretação e as fontes dos algoritmos estão em [docs/referencias-metodos.md](docs/referencias-metodos.md).

Um [painel inicial de proteínas](docs/proteinas.md) compara `recA`, `trhO` e `rpoB` nas quatro linhagens e mostra quais diferenças de códons preservam ou alteram os aminoácidos.

## Executar a análise

Os FASTA de `recA`, `trhO` e `rpoB` estão em `resources/<gene>/`. Eles podem ser recriados a partir dos genomas locais e das coordenadas em [proteins.csv](resources/proteins.csv):

```sh
python3 scripts/extract_proteins.py

make test
make

./build/gene-comparison --protein recA --benchmark
./build/gene-comparison --protein trhO --benchmark
./build/gene-comparison --protein rpoB --mode cds --benchmark
```

`--protein` é opcional e usa `recA` por padrão. O programa compara somente os seis pares únicos de cada gene em threads e constrói um alinhamento múltiplo progressivo das quatro proteínas. No modo `cds`, o alinhamento das proteínas orienta a projeção dos códons, preservando o quadro de leitura. `--threads N` limita a quantidade de threads; `--benchmark` também executa os pares sequencialmente e confere se os resultados são idênticos. Os arquivos gerados ficam em `results/<gene>/`: alinhamentos FASTA, tabela de pares, efeitos por códon e resumo de colunas. A [interpretação inicial](docs/resultados.md) reúne os três genes.
