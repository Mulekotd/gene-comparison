# Painel inicial de proteínas e diferenças de códons

Dados dos quatro genomas do projeto, com coordenadas e acessos RefSeq em [proteins.csv](../resources/proteins.csv). A presença abaixo significa que há uma **CDS anotada com produto proteico** no genoma; ela não demonstra expressão da proteína na célula. Os dados de anotação vieram dos arquivos GFF3 e FASTA de proteínas do [NCBI Genomes FTP](https://trace.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/). O projeto identifica a linhagem `1_2014MBL_MicDiv` como _M. dokdonensis_; o NCBI a nomeia _Maribacter_ sp.

| Proteína anotada                         | _M. aestuarii_   | _M. cobaltidurans_ | _M. dokdonensis_ (1_2014MBL_MicDiv) | _M. hydrothermalis_ | Presença |  Tamanho | Identidade entre proteínas |
| ---------------------------------------- | ---------------- | ------------------ | ----------------------------------- | ------------------- | -------: | -------: | -------------------------: |
| RecA — recombinase                       | `WP_281541862.1` | `WP_094996757.1`   | `WP_058105109.1`                    | `WP_068480290.1`    |      4/4 |   335 aa |               94,93–97,91% |
| TrhO — hidroxilase de uridina do tRNA    | `WP_281541863.1` | `WP_094996756.1`   | `WP_071340871.1`                    | `WP_068480292.1`    |      4/4 |   343 aa |               86,30–93,59% |
| RpoB — subunidade beta da RNA polimerase | `WP_281540825.1` | `WP_094997915.1`   | `WP_058106006.1`                    | `WP_068484266.1`    |      4/4 | 1.269 aa |               96,45–98,11% |

Cada gene aparece uma vez com esse nome na anotação de cada montagem. As quatro CDSs de cada linha traduzem exatamente para as proteínas RefSeq indicadas. Os seis alinhamentos entre pares de cada proteína não apresentaram lacunas. `recA` e `trhO` são vizinhos nos quatro genomas; `rpoB` mantém vizinhança com `rpoC`, `rplL`, `rplJ` e `rplA`. A [avaliação de ortologia](ortologia.md) reúne as evidências e seus limites para os três genes.

## Efeito observado das diferenças de bases

Comparações com `GCF_027474845.2` (_M. aestuarii_) como referência. As posições são alinhadas pelo códon e **excluem o códon de parada**. Uma diferença sinônima muda o códon e preserva o aminoácido; uma não sinônima muda o aminoácido codificado. A [explicação do código genético no NCBI](https://www.ncbi.nlm.nih.gov/sites/books/NBK21114/) descreve essas categorias.

| Gene   | Linhagem comparada  | Bases diferentes | Códons diferentes sem troca de aminoácido | Códons com troca de aminoácido |
| ------ | ------------------- | ---------------: | ----------------------------------------: | -----------------------------: |
| `recA` | 1_2014MBL_MicDiv    |              213 |                                       173 |                             13 |
| `recA` | _M. hydrothermalis_ |              199 |                                       167 |                             12 |
| `recA` | _M. cobaltidurans_  |              221 |                                       180 |                             17 |
| `trhO` | 1_2014MBL_MicDiv    |              245 |                                       162 |                             43 |
| `trhO` | _M. hydrothermalis_ |              232 |                                       148 |                             44 |
| `trhO` | _M. cobaltidurans_  |              249 |                                       160 |                             44 |
| `rpoB` | 1_2014MBL_MicDiv    |              596 |                                       488 |                             45 |
| `rpoB` | _M. hydrothermalis_ |              634 |                                       525 |                             45 |
| `rpoB` | _M. cobaltidurans_  |              619 |                                       527 |                             35 |

As três colunas numéricas usam unidades diferentes: **bases** na primeira e **posições de códons** nas outras duas. Por isso, não se somam diretamente; um único códon pode conter duas ou três bases diferentes. Como não há lacunas nesses alinhamentos, cada códon com troca de aminoácido corresponde a uma posição diferente na proteína.

| Exemplo, referência → 1_2014MBL_MicDiv | Posição na proteína | Códon         | Aminoácido | Interpretação          |
| -------------------------------------- | ------------------: | ------------- | ---------- | ---------------------- |
| `recA`                                 |                   2 | `AGC` → `AGT` | Ser → Ser  | Diferença sinônima     |
| `recA`                                 |                   6 | `GAC` → `GAA` | Asp → Glu  | Diferença não sinônima |
| `trhO`                                 |                  10 | `AAA` → `ATA` | Lys → Ile  | Diferença não sinônima |

Há também diferenças entre códons de parada, como `TAA` e `TAG`; ambos terminam a tradução e não entram na contagem de aminoácidos. As comparações descrevem **diferenças observadas entre linhagens**. Sem reconstruir a sequência ancestral, não é possível afirmar em qual linhagem cada alteração surgiu; uma troca de aminoácido, por si só, também não demonstra mudança de função.

Os números podem ser reproduzidos com `python3 scripts/analyze_proteins.py` na raiz do projeto. O script lê os genomas locais e as coordenadas de [proteins.csv](../resources/proteins.csv), extrai e traduz as CDSs e imprime as nove comparações em CSV.
