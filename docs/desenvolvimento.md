# Desenvolvimento

## Obtenção dos genomas

Todos os genomas utilizados no projeto foram obtidos na plataforma [NCBI](https://www.ncbi.nlm.nih.gov/).

## Próximas tarefas

- [x] **Preparar os dados:** selecionar os quatro genomas, documentar os acessos e as fontes e avaliar a [ortologia dos genes do painel](ortologia.md). A identificação de `GCF_001865655.1` como _M. dokdonensis_ foi confirmada para o projeto, embora o NCBI mostre _Maribacter_ sp. 1_2014MBL_MicDiv.
- [x] **Implementar o alinhamento e a análise de similaridade em C++:**
  - [x] Integrar as CDS e proteínas de `recA` em FASTA como entradas do programa; ler e validar as sequências, preservando o acesso e a linhagem de cada uma.
  - [x] Implementar o alinhamento global de duas sequências com Needleman–Wunsch, pontuação explícita e reconstrução do caminho; verificar o resultado com exemplos curtos de resposta conhecida.
  - [x] Gerar apenas os pares não ordenados de sequências (`i < j`). Para quatro linhagens, são `C(4,2) = 6` comparações; não repetir `B × A` após `A × B`.
  - [x] Processar os seis pares em paralelo com threads e uma fila de tarefas. Limitar o número de threads ao menor valor entre os núcleos disponíveis e os pares pendentes; manter cada resultado em posição própria para evitar disputas de escrita e apresentar a saída em ordem estável.
  - [x] Otimizar o alinhamento: trabalhar nas CDS e proteínas de `recA`, sem alinhar cromossomos completos; reutilizar sequências carregadas; calcular o escore com duas linhas da matriz e guardar apenas o necessário para reconstruir o alinhamento; evitar recalcular pares e métricas.
  - [x] Comparar os resultados da execução paralela com os da execução sequencial e medir o tempo de ambas. O programa aceita `--threads 1` para execução sequencial e `--benchmark` para comparar os tempos e conferir a igualdade dos resultados.
- [x] **Alinhar as quatro sequências:** o código constrói um alinhamento múltiplo progressivo das proteínas, preserva as colunas já introduzidas e projeta os códons das CDS sobre ele.
- [x] **Medir e interpretar:** os arquivos em `results/` contêm identidades entre pares, posições conservadas, substituições e lacunas; a [interpretação](resultados.md) compara DNA e proteína sem inferir função ou direção histórica apenas das diferenças observadas.
- [x] **Ampliar o programa para o [painel de proteínas](proteinas.md):** `--protein` seleciona `recA`, `trhO` ou `rpoB`; o programa calcula seis pares únicos e apresenta diferenças sinônimas e não sinônimas por códon. As contagens foram comparadas com o script independente do painel.
- [ ] **Entregar:** publicar o código, os dados de entrada com seus acessos, os alinhamentos, as tabelas de resultados e o artigo no modelo SBC.

**Implementação concluída para o painel inicial.** Falta a entrega final do artigo no modelo SBC e a revisão dos resultados biológicos nele apresentados.
