# Avaliação de ortologia de `recA`

Consulta realizada em 24 de novembro de 2026. A ortologia se aplica aos genes `recA`, não aos genomas inteiros.

## Dados examinados

Foram usados os quatro genomas FASTA do projeto e, para cada montagem, os arquivos RefSeq `genomic.gff.gz` e `protein.faa.gz` do [NCBI Genomes FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/). As posições abaixo são inclusivas, em coordenadas de 1 a N.

| Montagem | Linhagem | Cromossomo | CDS `recA` | Fita | Locus tag | Proteína |
| --- | --- | --- | ---: | :---: | --- | --- |
| `GCF_027474845.2` | _M. aestuarii_ JCM18631 | `NZ_CP107031.2` | 1.855.852–1.856.859 | − | `N8A89_RS08400` | [`WP_281541862.1`](https://www.ncbi.nlm.nih.gov/protein/WP_281541862.1) |
| `GCF_002269385.1` | _M. cobaltidurans_ B1 | `NZ_CP022957.1` | 1.787.383–1.788.390 | + | `CJ263_RS07810` | [`WP_094996757.1`](https://www.ncbi.nlm.nih.gov/protein/WP_094996757.1) |
| `GCF_001865655.1` | _Maribacter_ sp. 1_2014MBL_MicDiv | `NZ_CP011318.1` | 2.060.344–2.061.351 | − | `YQ22_RS08670` | [`WP_058105109.1`](https://www.ncbi.nlm.nih.gov/protein/WP_058105109.1) |
| `GCF_001913155.1` | _M. hydrothermalis_ T28 | `NZ_CP018760.1` | 2.058.793–2.059.800 | − | `BTR34_RS08805` | [`WP_068480290.1`](https://www.ncbi.nlm.nih.gov/protein/WP_068480290.1) |

**Nomenclatura:** o arquivo `samples.csv` classifica `GCF_001865655.1` como _M. dokdonensis_, enquanto o cabeçalho do [registro cromossômico `NZ_CP011318.1`](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP011318.1) e a anotação RefSeq usam _Maribacter_ sp. 1_2014MBL_MicDiv. O responsável pelo projeto confirmou a identificação como _M. dokdonensis_; convém registrar essa diferença de nomes no relatório final.

## Verificações

- Cada anotação GFF contém **uma** CDS com `gene=recA`. As quatro CDSs têm 1.008 nucleotídeos, começam com `ATG`, terminam em códon de parada e traduzem exatamente para as proteínas RefSeq correspondentes, de 335 aminoácidos, sem parada interna.
- As seis comparações globais entre proteínas têm 335 posições alinhadas, **sem lacunas**. Identidades: aestuarii–cobaltidurans **94,93%** (318/335); aestuarii–linhagem 1_2014MBL_MicDiv **96,12%** (322/335); aestuarii–hydrothermalis **96,42%** (323/335); cobaltidurans–linhagem 1_2014MBL_MicDiv **95,52%** (320/335); cobaltidurans–hydrothermalis **95,52%** (320/335); linhagem 1_2014MBL_MicDiv–hydrothermalis **97,91%** (328/335).
- Uma triagem de todas as proteínas anotadas em cada montagem, usando sobreposição de segmentos exatos de quatro aminoácidos, colocou a RecA esperada em primeiro lugar nas **12 comparações direcionais** entre proteomas. Ela compartilhou de 274 a 303 segmentos com a RecA consultada; a segunda colocada compartilhou no máximo 25. Essa triagem é uma verificação independente da anotação, não uma busca BLAST formal.
- A vizinhança gênica é conservada nas quatro montagens. Considerando o sentido de transcrição de `recA`, aparecem `trhO`, uma proteína com domínio PSP1, `gldH` e uma proteína ligadora de penicilina 1A de um lado; uma proteína hipotética e um fator sigma da RNA polimerase aparecem do outro. A orientação do bloco no cromossomo varia entre montagens.

## Conclusão

As evidências sustentam **fortemente** que as quatro CDSs de `recA` são ortólogas: uma cópia anotada por genoma, proteínas completas e muito semelhantes, correspondência destacada entre proteomas e conservação da vizinhança. A classificação evolutiva não fica demonstrada de modo absoluto por esses critérios; uma análise filogenética com linhagens adicionais poderia investigar duplicações antigas ou transferência horizontal.

As identidades foram calculadas por alinhamento global com pontuação +2 para coincidência, −1 para substituição e −3 para lacuna. A distinção entre ortologia e paralogia segue a [definição do NCBI](https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/ortho_para.html); limitações de correspondência recíproca em bactérias são discutidas em [Bacterial Microevolution and the Pangenome](https://www.ncbi.nlm.nih.gov/books/NBK558819/).
