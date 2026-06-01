# CellDeath_Figure 已选图件指标与坐标轴说明

本文件说明 `CellDeath_Figure` 中已选主图和补充图的绘图指标、横纵坐标含义和必要的解释边界。所有计数均为文献层面的 PubMed/PubTator 记录、共现记录或由这些记录派生的统计量。请将这些图理解为文献图谱、术语使用、证据阶段或文献结构信号；它们不等同于因果生物学证明、临床疗效证据或通路实验验证。

## 通用指标

- `PMID`：PubMed 文献唯一编号。唯一 PMID 数表示相关过滤步骤后保留的不同 PubMed 记录数。
- `PMID-stage records`：与证据阶段、疾病、死亡方式、基因或药物语境相连的 PMID 层级记录。一篇文献如果支持多个语境，可以出现在多行记录中。
- `Pair_Count`：死亡方式 x 疾病共现对的文献支持计数。这里是文献支持量，不是实验样本量。
- `Disease_Breadth`：过滤后与某一死亡方式相连的疾病术语数量。
- `Death_Mode`：图谱使用的调控性细胞死亡概念标签。本组选图围绕 8 类最终入选的死亡概念：ferroptosis、pyroptosis、necroptosis、NETosis、immunogenic cell death、cuproptosis、PANoptosis 和 disulfidptosis。
- `log2 RCA` / `log2 selectivity`：log2 相对集中度或选择性指标。正值表示相对背景富集，负值表示低于背景。
- `hub score`：由文献体量、疾病广度、PMID 支持和相关文献结构特征构成的复合文献枢纽分数。该指标用于排序文献枢纽，不代表生物学枢纽机制。
- `multi_death_signal_score`：基因层面的多死亡方式文献信号分数，定义为 `death_mode_breadth x log10(PMID count + 1)`。
- `oncology fraction`：某一基因、药物或死亡方式文献信号中，与肿瘤/肿瘤语境相关的比例。
- `approved-drug` 层：从 PubMed/PubTator 药物实体出发，映射到已获批临床药物参考表后的药物文献信号。绘图前排除了常见离子、营养物质和代谢物样术语。

## Figure 1

**Panel A. 疾病文献中的年度发文量**
- 横坐标：`Year`，2000-2025。
- 纵坐标：年度文献数量，即每个死亡方式在疾病文献中的唯一 PMID 数。
- 颜色/折线：死亡方式。
- 背景阴影：用于叙事解释的术语采用阶段。

**Panel B. 2000-2025 过滤后的起始年、起飞年和峰值年**
- 横坐标：年份。
- 纵坐标：死亡方式。
- 空心圆：该死亡方式在过滤后图谱中的首次出现年份。
- 彩色圆：起飞年，由时间趋势汇总逻辑定义。
- 黑色方块：峰值年，即年度发文量最高的年份。

主要源表：`table/01_Figure_1_temporal_adoption/death_yearly_selected_death_modes_2000_2025.csv` 和 `temporal_summary_selected_death_modes_2000_2025.csv`。

## Figure 2

**Panel A. 过滤后的文献广度与文献体量不均衡**
- 横坐标：死亡方式唯一 PMID 数，通常为 log 尺度。
- 纵坐标：疾病广度，即与该死亡方式相连的疾病术语数量。
- 点大小/颜色：死亡方式的文献体量和类别标识。

**Panel B. 原始疾病枢纽**
- 横坐标：文献枢纽分数。
- 纵坐标：去除肿瘤前的高排名疾病枢纽术语。
- 条形颜色：疾病系统类别或肿瘤/非肿瘤类别。
- 标签：用于计算枢纽分数的疾病广度和 pair/PMID 支持信息。

**Panel C. 国家/地区层面的署名单位信号**
- 横坐标：PMID 数。
- 纵坐标：国家/地区或署名单位标签。
- 堆叠颜色：该国家/地区信号中不同死亡方式的贡献。
- 解释边界：表示署名单位层面的发文语境，不代表作者国籍，也不是国家科研质量排名。

**Panel D. 图谱文献数最高的期刊**
- 横坐标：PMID 数。
- 纵坐标：期刊。
- 堆叠颜色：每个期刊中不同死亡方式的贡献。
- 解释边界：表示图谱文献在期刊中的集中度，不是期刊质量指标。

**Panel E. 疾病语境雷达图**
- 每个小雷达图：一个死亡方式。
- 角向标签：该死亡方式内选出的高排名疾病术语。
- 半径：来自 balanced top-pair 源表的证据阶段信号。
- 解释边界：用于压缩展示疾病语境画像，不代表临床成熟度证明。

主要源目录：`table/02_Figure_2_breadth_hubs_journal_country_radar/`。

## Figure 3

**Panel A. 审计后的证据阶段构成**
- 横坐标：PMID-stage records 的比例，范围 0-1。
- 纵坐标：死亡方式。
- 堆叠填充：审计后的证据阶段类别，包括 `Basic mechanism`、`Preclinical animal`、`Human observational`、`PublicationType-supported trial` 以及 guideline-like 不合格记录。
- 解释边界：分母是 PMID-stage records，不是唯一研究数。

**Panel B. 审计后的高阶段证据信号**
- 横坐标：比例。
- 纵坐标：死亡方式。
- 灰色条形：原始文本高阶段信号。
- 蓝色条形/点：经过 PublicationType 支持后的 trial 信号。
- 注释：审计计数和降级规则，包括 soft clinical terms 的重新分类。

主要源目录：`table/03_Figure_3_evidence_stage_audit/`。

## Figure 4

**Panel A. Top 40 广谱多死亡方式基因信号**
- 横坐标：支持该共享基因的 PMID 数。
- 纵坐标：共享基因符号。
- 填充/文字：死亡方式广度，即该基因连接的死亡方式数量。
- 排序：`multi_death_signal_score = death_mode_breadth x log10(PMID count + 1)`。

**Panel B. 基因 x 死亡方式热图**
- 横坐标：死亡方式。
- 纵坐标：基因。
- 填充：该基因 x 死亡方式组合的 `log10(PMID count + 1)` 文献支持强度。

**Panel C. 肿瘤富集基因**
- 横坐标：PMID 数。
- 纵坐标：基因。
- 填充：多死亡方式信号分数。
- 纳入条件：具有较高 oncology fraction 的基因。

**Panel D. 非肿瘤富集基因**
- 横坐标：PMID 数。
- 纵坐标：基因。
- 填充：多死亡方式信号分数。
- 纳入条件：具有较高 non-oncology fraction 的基因。

主要源目录：`table/04_Figure_4_gene_literature_layer/`。

## Figure 5

**Panel A. 肿瘤文献热度指数**
- 横坐标：分量整合后的 heat index。
- 纵坐标：死亡方式。
- 条形颜色：死亡方式。
- heat index 组成：肿瘤唯一 PMIDs、肿瘤家族广度、近期增长，以及 refined approved-drug PMID fraction。

**Panel B. 肿瘤比例与已获批药物词汇信号**
- 横坐标：死亡方式记录中属于肿瘤语境的比例。
- 纵坐标：approved-drug PMID fraction。
- 点大小：refined approved clinical drug entities 数量。
- 点标签/颜色：死亡方式。

**Panel C. 肿瘤家族 x 死亡方式选择性**
- 横坐标：死亡方式。
- 纵坐标：肿瘤家族。
- 填充：经过截断的 `log2_tumor_death_selectivity`，相对于肿瘤背景归一化。
- 点大小：pair count。

**Panel D. 参考药物类别/治疗类别 x 死亡方式候选景观**
- 横坐标：死亡方式。
- 纵坐标：参考类别或治疗类别。
- 填充：经过截断的 log2 selectivity。
- 点大小：候选药物或关系支持量，具体取决于源表。

**Panel E. 最高支持的肿瘤-药物-死亡方式三元组**
- 横坐标：唯一 PMID 数。
- 纵坐标：肿瘤家族 | 药物 | 死亡方式三元组。
- 条形填充：死亡方式。
- 白点：induction/sensitization 句子候选 PMID 支持量。

主要源目录：`table/05_Figure_5_oncology_therapeutic_landscape/`。

## Figure S1

**疾病文献中的归一化年度轨迹**
- 横坐标：年份。
- 纵坐标：归一化年度计数；每个死亡方式内部缩放，使其最大年度值等于 1。
- 颜色/折线：死亡方式。
- 背景阴影：术语采用阶段。

主要源目录：`table/06_Figure_S1_normalized_temporal/`。

## Figure S2

**Panel A. 低计数关系仍然数量较多**
- 横坐标：pair count，绘图时在 50 处截断。
- 纵坐标：pair 数量。
- 垂直参考线：稳健性分析中使用的 pair-count 阈值。
- 注：最后一个 bin 合并所有大于或等于截断值的 pair。

**Panel B. 阈值保留情况**
- 横坐标：`Pair_Count` 阈值。
- 纵坐标：保留数量。
- 折线颜色：保留 pair 或保留 disease 的数量。

主要源目录：`table/07_Figure_S2_low_count_threshold/`。

## Figure S3

**Panel A. 疾病系统构成**
- 横坐标：positive pairs 数量。
- 纵坐标：疾病系统。
- 条形颜色：疾病系统类别。

**Panel B. 去除肿瘤后的疾病枢纽**
- 横坐标：文献枢纽分数。
- 纵坐标：保留的非肿瘤疾病枢纽。
- 解释边界：用于评估肿瘤文献体量敏感性，不是生物学枢纽证明。

**Panel C. 图谱文献数最高的期刊**
- 横坐标：图谱文献中的死亡方式比例或计数，取决于具体面板源表。
- 纵坐标：期刊。
- 堆叠填充：死亡方式。

**Panel D. 国家/地区层面的署名单位信号**
- 横坐标：图谱文献中各死亡方式的比例。
- 纵坐标：国家/地区或署名单位标签。
- 堆叠填充：死亡方式。

**Panel E. 高影响期刊中的国家/地区署名单位信号**
- 横坐标：通过 high-impact/JIF 匹配规则后的文章级 PMID 数。
- 纵坐标：国家/地区或署名单位标签。
- 堆叠填充：死亡方式。

**Panel F. 分母归一化后的期刊语境**
- 横坐标：journal atlas share 或 RCD-focus index，由 PubMed 分母归一化得到。
- 纵坐标：期刊。
- 解释边界：表示发表生态语境，不是期刊质量评价。

主要源目录：`table/08_Figure_S3_ecosystem_and_hubs/`。

## Figure S4

**Panel A. 国家/地区专门化热图**
- 横坐标：死亡方式。
- 纵坐标：国家/地区。
- 填充：log2 RCA / specialization score。
- 正值含义：该国家/地区信号相对于图谱背景更富集于对应死亡方式。

**Panel B. 期刊专门化热图**
- 横坐标：死亡方式。
- 纵坐标：期刊。
- 填充：log2 RCA / specialization score。
- 正值含义：该期刊信号相对于图谱背景更富集于对应死亡方式。

主要源目录：`table/09_Figure_S4_country_journal_heatmaps/`。

## Figure S5

**Panel A. 标题/摘要/关键词贡献**
- 横坐标：匹配字段，即 `title`、`abstract` 或 `keywords`。
- 纵坐标：primary gene co-mention matches 数量。

**Panel B. Top 死亡方式选择性基因信号**
- 横坐标：PMID 数。
- 纵坐标：基因。
- 填充：该基因具有选择性的死亡方式。
- 纳入条件：死亡方式广度等于 1。

**Panel C. STRING/network 面板**
- 网络节点：STRING 导出中使用的共享基因或蛋白。
- 边：STRING 蛋白互作或导出的网络连接。
- 节点/边样式：来自外部 STRING/VOSviewer 风格导出；已复制的 SVG/PNG 是该面板的视觉源资产。

**Panel D. 肿瘤比例与多死亡方式分数**
- 横坐标：oncology fraction。
- 纵坐标：multi-death signal score。
- 颜色：富集类别，包括 `oncology_enriched`、`non_oncology_enriched` 或 `mixed`。

主要源目录：`table/10_Figure_S5_gene_workflow_string_scatter/`。

## Figure S6

**Panel A. 非肿瘤 VOSviewer 疾病-基因网络**
- 节点位置：VOSviewer 布局坐标 `x` 和 `y`。
- 节点标签：疾病或基因条目标签。
- 节点聚类/颜色：VOSviewer 或 Louvain 聚类编号。
- 节点权重字段：`weight<Genes>`、`weight<PMIDs>` 和 `weight<DeathModes>`。
- 边权重：网络文件中由 PMID 支持的共现权重。
- 过滤：非肿瘤图谱按 `weight<PMIDs> > 30` 过滤并重新聚类。

**Panel B. 肿瘤 VOSviewer 疾病-基因网络**
- 坐标、聚类、节点权重和边权重含义与 Panel A 相同。
- 范围：肿瘤疾病-基因文献网络。

主要源目录：`table/11_Figure_S6_vosviewer_networks/`。

## Figure S7

**Panel A. Top 已获批临床药物 x 死亡方式**
- 横坐标：死亡方式。
- 纵坐标：已获批临床药物。
- 填充：`log1p(PMID count)` 或 PMID 支持强度。
- 点大小：induction/sensitization candidate PMID 支持量。

**Panel B. 各死亡方式的已获批药物 PMID 支持**
- 横坐标：汇总后的 drug-mode PMID 支持量。
- 纵坐标：死亡方式。
- 条形颜色：死亡方式。

**Panel C. 已获批临床药物按死亡方式分层的文献支持**
- 横坐标：唯一 PMID 数。
- 纵坐标：已获批临床药物。
- 填充/分面：死亡方式。
- 解释边界：表示肿瘤死亡方式语境中已获批药物提及的文献支持，不代表疗效证据。

主要源目录：`table/12_Figure_S7_approved_drug_landscape/`。
