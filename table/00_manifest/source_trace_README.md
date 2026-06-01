# CellDeath_Figure 源表、代码与 raw PMID 溯源说明

本目录汇总了 `CellDeath_Figure` 已选图件从上游结果到绘图所需的源表、脚本和最原始 PMID/PubMed/PubTator 层文件，便于后续复现、审稿补充和人工核查。

当前整理结果：

- `table`：112 个文件，按图号组织源表、绘图中间结果和初筛汇总表。
- `code`：44 个文件，按分析分支组织绘图和上游处理脚本。
- `raw`：9 个 PMID/PubMed/PubTator 语境文件。
- 已从整理后的 `table/`、`code/` 和 `raw/` 目录中移除 macOS `._*` sidecar 文件。

## Table 目录

- `table/00_manifest/selected_death_modes_manifest.csv`：最终 8 种细胞死亡方式的 PubMed term-screen 文章数、进入疾病 pair 图谱的唯一 PMID 数和 pair-atlas 中首次出现年份。
- `table/01_Figure_1_temporal_adoption`：年度 PMID 计数和起飞年/峰值年时间趋势汇总。
- `table/02_Figure_2_breadth_hubs_journal_country_radar`：文献广度-体量、疾病枢纽、期刊/国家、高影响期刊国家信号和雷达图源表。
- `table/03_Figure_3_evidence_stage_audit`：证据阶段审计、soft clinical terms 降级和高阶段信号源表。
- `table/04_Figure_4_gene_literature_layer`：gene-v3 广谱基因、选择性基因、肿瘤/非肿瘤基因层源表。
- `table/05_Figure_5_oncology_therapeutic_landscape`：肿瘤热度、肿瘤比例、肿瘤家族选择性、治疗类别和肿瘤-药物-死亡方式三元组源表。
- `table/06_Figure_S1_normalized_temporal`：归一化年度轨迹源表。
- `table/07_Figure_S2_low_count_threshold`：低计数关系和阈值保留稳健性源表。
- `table/08_Figure_S3_ecosystem_and_hubs`：疾病系统、去肿瘤枢纽、期刊、国家和分母归一化发表生态源表。
- `table/09_Figure_S4_country_journal_heatmaps`：国家/期刊 x 死亡方式计数矩阵和 log2 RCA 矩阵。
- `table/10_Figure_S5_gene_workflow_string_scatter`：基因共现流程、STRING 视觉资产、选择性基因和肿瘤比例散点图源表。
- `table/11_Figure_S6_vosviewer_networks`：VOSviewer map/network/count/source-summary 文件，以及非肿瘤和肿瘤网络图源资产。
- `table/12_Figure_S7_approved_drug_landscape`：已获批临床药物景观、参考药物表、药物映射和关系句子源表。

## Code 目录

- `code/01_selected_death_modes_v5`：v5 selected death-mode 文献图谱主线的 Python 绘图和整理脚本。
- `code/02_gene_v3`：gene-v3 扩展分支的 R 脚本，以及 `comorbidity_pubmed_pipeline.py`。
- `code/03_oncology_therapeutic_landscape`：肿瘤治疗/药物景观分支的 R 与 Python 脚本。

## Raw PMID 目录

- `raw/01_core_pubmed_atlas`：核心死亡方式-疾病 pair-PMID 文件，以及 selected death-mode PMID-stage atlas records。
- `raw/02_gene_v3_pubmed_context`：gene-v3 PMID 文本语境、pair-PMID links 和原始 gene co-mention matches。
- `raw/03_oncology_pubmed_pubtator`：肿瘤文章记录、PubMed/PubTator 药物语境、候选药物映射和未过滤候选语境。

## 图件来源范围

- `Figure 1`、`Figure 2`、`Figure 3`、`Figure S1`、`Figure S2`、`Figure S3` 和 `Figure S4` 主要来自 v5 selected death-mode 文献图谱分支。
- `Figure 4`、`Figure S5` 和 `Figure S6` 来自 gene-v3 / VOSviewer 扩展分支。
- `Figure 5` 和 `Figure S7` 来自 oncology therapeutic / approved-drug landscape 分支。

## 使用边界

- `table` 中的文件用于追踪从上游统计结果到最终绘图的数据输入，不代表所有原始全文或数据库快照。
- `raw` 中保留的是最接近上游的 PMID/PubMed/PubTator 语境文件，可用于追溯每个图件背后的文献记录。
- 所有图件指标均为文献层级信号，应按文献图谱、证据阶段或语境富集解释，不应写成因果机制或临床疗效结论。
