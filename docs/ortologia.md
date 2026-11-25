# Avaliação de ortologia dos genes do painel

Esta avaliação examina `recA`, `trhO` e `rpoB` nas quatro linhagens de _Maribacter_ do projeto. **Ortologia é uma relação entre genes de linhagens diferentes**, ligada à sua história evolutiva; semelhança de sequência ou igualdade do nome anotado, isoladamente, não a demonstram. A distinção entre ortólogos e parálogos é explicada pelo [NCBI](https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/ortho_para.html) e por [Tatusov, Koonin e Lipman (1997)](https://pubmed.ncbi.nlm.nih.gov/9381173/).

Foram usados os genomas FASTA locais e, para cada montagem, os arquivos RefSeq `genomic.gff.gz` e `protein.faa.gz` do [NCBI Genomes FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/). Os acessos das proteínas, cromossomos e coordenadas inclusivas das 12 CDSs estão em [proteins.csv](../resources/proteins.csv). As sequências extraídas estão em `resources/<gene>/`.

## Evidências nos três genes

| Gene   | Cópias anotadas | CDS por linhagem | Proteína por linhagem | Identidade entre os seis pares | Primeira correspondência na triagem entre proteomas | Vizinhança conservada                                                                       |
| ------ | --------------: | ---------------: | --------------------: | -----------------------------: | --------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| `recA` |  1 por montagem |      1.008 bases |                335 aa |                   94,93–97,91% | 12 de 12 buscas direcionais                         | `trhO`, proteína com domínio PSP1, `gldH`, proteína ligadora de penicilina 1A e fator sigma |
| `trhO` |  1 por montagem |      1.032 bases |                343 aa |                   86,30–93,59% | 12 de 12 buscas direcionais                         | `recA`, proteína com domínio PSP1 e `gldH`                                                  |
| `rpoB` |  1 por montagem |      3.810 bases |              1.269 aa |                   96,45–98,11% | 12 de 12 buscas direcionais                         | `rpoC`, `rplL`, `rplJ` e `rplA`                                                             |

As 12 CDSs começam com `ATG`, terminam em códon de parada, não têm parada interna e traduzem exatamente para os produtos RefSeq indicados em `proteins.csv`. Os alinhamentos globais das proteínas não apresentaram lacunas; os resultados por par estão em `results/<gene>/pairs.tsv`.

A triagem comparou cada proteína do painel com **todas as proteínas anotadas** nos outros três proteomas, usando o número de segmentos distintos e idênticos de quatro aminoácidos. A proteína esperada ficou em primeiro lugar em todas as 12 buscas direcionais de cada gene. As faixas de segmentos compartilhados com a primeira colocada foram 274–303 (`recA`), 200–270 (`trhO`) e 1.107–1.179 (`rpoB`); a segunda colocada compartilhou, no máximo, 25, 20 e 82 segmentos, respectivamente. Trata-se de uma triagem heurística, não de BLAST recíproco formal nem de prova filogenética.

A ordem dos genes vizinhos também se conserva. Em `GCF_002269385.1`, os blocos aparecem na orientação cromossômica oposta à de algumas outras montagens; isso não altera a correspondência dos genes.

## Exemplo detalhado: `recA`

| Montagem          | Linhagem                | Cromossomo      |          CDS `recA` | Fita | Locus tag       | Proteína                                                                |
| ----------------- | ----------------------- | --------------- | ------------------: | :--: | --------------- | ----------------------------------------------------------------------- |
| `GCF_027474845.2` | _M. aestuarii_ JCM18631 | `NZ_CP107031.2` | 1.855.852–1.856.859 |  −   | `N8A89_RS08400` | [`WP_281541862.1`](https://www.ncbi.nlm.nih.gov/protein/WP_281541862.1) |
| `GCF_002269385.1` | _M. cobaltidurans_ B1   | `NZ_CP022957.1` | 1.787.383–1.788.390 |  +   | `CJ263_RS07810` | [`WP_094996757.1`](https://www.ncbi.nlm.nih.gov/protein/WP_094996757.1) |
| `GCF_001865655.1` | 1_2014MBL_MicDiv        | `NZ_CP011318.1` | 2.060.344–2.061.351 |  −   | `YQ22_RS08670`  | [`WP_058105109.1`](https://www.ncbi.nlm.nih.gov/protein/WP_058105109.1) |
| `GCF_001913155.1` | _M. hydrothermalis_ T28 | `NZ_CP018760.1` | 2.058.793–2.059.800 |  −   | `BTR34_RS08805` | [`WP_068480290.1`](https://www.ncbi.nlm.nih.gov/protein/WP_068480290.1) |

As seis identidades de RecA são: aestuarii–cobaltidurans **94,93%** (318/335), aestuarii–1_2014MBL_MicDiv **96,12%** (322/335), aestuarii–hydrothermalis **96,42%** (323/335), cobaltidurans–1_2014MBL_MicDiv **95,52%** (320/335), cobaltidurans–hydrothermalis **95,52%** (320/335) e 1_2014MBL_MicDiv–hydrothermalis **97,91%** (328/335).

**Nomenclatura:** a planilha do projeto classifica `GCF_001865655.1` como _M. dokdonensis_, enquanto o [registro cromossômico `NZ_CP011318.1`](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP011318.1) e a anotação RefSeq usam _Maribacter_ sp. 1*2014MBL_MicDiv. O responsável pelo projeto confirmou a identificação como \_M. dokdonensis*; essa diferença de nomes deve constar no relatório final.

## Conclusão e limites

A cópia única anotada, a tradução íntegra, a alta semelhança das proteínas, a correspondência destacada entre proteomas e a vizinhança conservada **sustentam fortemente a hipótese de ortologia** para os três genes nas quatro linhagens. Isso ainda é uma inferência: duplicações antigas seguidas de perda de cópias ou transferência horizontal podem confundir a relação evolutiva. Uma árvore de cada gene com linhagens adicionais, comparada a uma árvore de espécies, fortaleceria a avaliação. O NCBI também descreve os [limites de buscas de melhores correspondências em bactérias](https://www.ncbi.nlm.nih.gov/books/NBK558819/).

As identidades foram calculadas por alinhamento global com pontuação `+2` para coincidência, `−1` para substituição e `−3` para lacuna. A [interpretação das diferenças](resultados.md) não atribui direção histórica a substituições ou lacunas. As [referências de método](referencias-metodos.md) separam a evidência biológica das escolhas do programa.
