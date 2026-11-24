# Desenvolvimento

## Obtenção dos genomas

Todos os genomas utilizados no projeto foram obtidos na plataforma [NCBI](https://www.ncbi.nlm.nih.gov/).

## Próximas tarefas

- [x] **Preparar os dados:** selecionar os quatro genomas, documentar os acessos e as fontes e avaliar a [ortologia de `recA`](ortologia-reca.md). A identificação de `GCF_001865655.1` como _M. dokdonensis_ foi confirmada para o projeto, embora o NCBI mostre _Maribacter_ sp. 1_2014MBL_MicDiv.
- [ ] **Implementar o alinhamento e a análise de similaridade em C++:**
  - [ ] Integrar as CDS e proteínas de `recA` em FASTA como entradas do programa; ler e validar as sequências, preservando o acesso e a linhagem de cada uma.
  - [ ] Implementar o alinhamento global de duas sequências com Needleman–Wunsch, pontuação explícita e reconstrução do caminho; verificar o resultado com exemplos curtos de resposta conhecida.
  - [ ] Gerar apenas os pares não ordenados de sequências (`i < j`). Para quatro linhagens, são `C(4,2) = 6` comparações; não repetir `B × A` após `A × B`.
  - [ ] Processar os seis pares em paralelo com threads e uma fila de tarefas. Limitar o número de threads ao menor valor entre os núcleos disponíveis e os pares pendentes; manter cada resultado em posição própria para evitar disputas de escrita e apresentar a saída em ordem estável.
  - [ ] Otimizar o alinhamento: trabalhar nas CDS e proteínas de `recA`, sem alinhar cromossomos completos; reutilizar sequências carregadas; calcular o escore com duas linhas da matriz e guardar apenas o necessário para reconstruir o alinhamento; evitar recalcular pares e métricas.
  - [ ] Comparar os resultados da execução paralela com os da execução sequencial e medir o tempo de ambas. Usar a versão paralela quando houver ganho real para o volume de dados.
- [ ] **Alinhar as quatro sequências:** estender o código para um alinhamento múltiplo progressivo que preserve as colunas e lacunas já introduzidas. Alinhar as proteínas e projetar o resultado sobre as CDS para comparar códons.
- [ ] **Medir e interpretar:** calcular identidades entre pares, posições conservadas, substituições e lacunas; comparar as alterações no DNA e na proteína. Evitar inferir a direção histórica de uma lacuna isolada ou mudança de função apenas por uma substituição de aminoácido.
- [ ] **Entregar:** publicar o código, os dados de entrada com seus acessos, os alinhamentos, as tabelas de resultados e o artigo no modelo SBC.

**Próximo marco:** ter um alinhador global de duas sequências funcionando com exemplos pequenos e as seis comparações únicas reproduzíveis em execução sequencial e paralela. Relatar os resultados biológicos após essa validação.
