# ABCB5 SNP Analysis — Upstream Pipeline (GATK)

> 一套**可复现、可嵌入任意数据**的 germline SNP 上游检测与分型流程（基于 **GATK HaplotypeCaller**）。
> 从原始测序数据（FASTQ/BAM）到高质量 SNP 位点集、基因型矩阵与 lollipop 可视化，
> 样本、参考、区域与阈值均由配置文件决定，不锁定任何具体样本或结果。
>
> **本仓库仅包含上游 pipeline**（建索引 → 比对 → HaplotypeCaller → 过滤 → 合并 → 基因型矩阵 → 可视化）。
> 下游分析（双等位金标准提纯、下采样召回率、与 bcftools 对比等）不在此仓库中。
> bcftools 版本见姊妹仓库 **ABCB5-SNP**。

---

## 目录

- [概览](#概览)
- [分析流程](#分析流程)
- [样本](#样本)
- [文件结构](#文件结构)
- [快速使用](#快速使用)
- [依赖环境](#依赖环境)
- [结果说明（含"表格错位"修复说明）](#结果说明含表格错位修复说明)

---

## 概览

本仓库交付一套**完整、可复现**的上游分析 pipeline，用于从原始测序数据
（FASTQ 或已比对 BAM）中检出 ABCB5 区域内的 germline SNP，并产出：

- 过滤后的高质量双等位 SNP 位点集（合并多样本 VCF）
- 各样本基因型矩阵
- 各样本 SNP 位点分布 lollipop 图 + 位点对照表（`variant_positions.csv`，Excel 可直接分列）

## 分析流程

| 步骤 | 工具 | 说明 |
|------|------|------|
| 建索引 | `bwa index` / `samtools faidx` / `gatk CreateSequenceDictionary` | 参考基因组建索引（一次性） |
| 比对 | BWA-MEM | 基于用户提供的参考基因组与目标区域 BED（可选，已有 BAM 可跳过） |
| 变异检测 | `gatk HaplotypeCaller` | 逐样本 calling，ploidy 可配 |
| 过滤 | `bcftools view` | 双等位 SNP，DP≥30, MQ≥30, VAF≥0.30, GT 为 0/1 或 1/1 |
| 合并 | `bcftools merge` | 多样本合并为统一 VCF |
| 基因型矩阵 | `bcftools query` | 输出 `%CHROM %POS %REF %ALT %GT` 矩阵 |
| 可视化 | `plot_lollipop_gatk.py` | 各样本 SNP 分布图 + `variant_positions.csv`（可选，需 domains 文件） |

## 样本

样本由 `config/sample_info.csv` 决定，**不写死在代码中**。
在该 CSV 中按 `sample_name,group` 格式逐行列出你的样本即可；
样本名须与 `input/` 目录下 `<样本名>_R1.fq.gz` / `<样本名>_R2.fq.gz` 对应（或已有 BAM 时与 `alignment/<样本>.sorted.bam` 对应）。
（仓库自带一份示例样本清单，仅作模板，可替换为任意数据。）

若你已用其他方式完成比对，可设 `SKIP_ALIGN=true` 跳过比对步骤，
pipeline 将直接读取 `results/alignment/<样本>.sorted.bam`。

## 文件结构

```
ABCB5-SNP-GATK/
├── README.md
├── .gitignore
├── scripts/
│   ├── config.sh              # 所有路径/阈值外部化配置
│   ├── run_pipeline.sh         # 上游 pipeline 主脚本（一键复现）
│   └── plot_lollipop_gatk.py   # lollipop 图 + 位点对照表（CSV）
├── config/
│   ├── sample_info.csv         # 样本清单（用户自定义）
│   └── domains.txt             # 蛋白结构域坐标（lollipop 图用，可选）
├── data/                       # 用户自备：reference.fa / targets.bed
├── input/                      # 用户自备：<样本名>_R1.fq.gz / _R2.fq.gz
└── results/                    # 运行产物（不纳入版本控制）
    ├── alignment/              # 比对 BAM
    └── gatk/upstream/          # 各样本 raw/filtered VCF + combined + genotypes + lollipop/
```

> 运行产物（`results/`、`data/`、`input/`）由 pipeline 动态生成或由用户提供，
> 已列入 `.gitignore`，不会随仓库分发。

## 快速使用

```bash
# 1. 克隆仓库
git clone https://github.com/18823534222-dev/ABCB5-SNP-GATK.git
cd ABCB5-SNP-GATK

# 2. 准备用户输入
#    - data/reference.fa, data/targets.bed
#    - input/<样本名>_R1.fq.gz / _R2.fq.gz
#    - 编辑 config/sample_info.csv 列出你的样本
#    - （可选）编辑 scripts/config.sh 调整 GATK 路径与阈值

# 3. 激活 conda 环境（含 bwa / samtools / bcftools / gatk / python+matplotlib）
conda activate snpcaller

# 4. 运行上游 pipeline
bash scripts/run_pipeline.sh

# 已有比对好的 BAM 时，跳过比对：
SKIP_ALIGN=true bash scripts/run_pipeline.sh
```

> 注：仓库本身**不包含运行产物与用户数据**（已列入 `.gitignore`）。
> 克隆后请自备 `data/`（参考与 BED）与 `input/`（样本 FASTQ，或已比对 BAM），
> 运行 `bash scripts/run_pipeline.sh` 即可在本机生成 `results/` 下全部结果。

## 依赖环境

- **比对**：`bwa` (BWA-MEM)
- **处理**：`samtools`, `bedtools`
- **变异检测**：`gatk` (HaplotypeCaller, 4.x)
- **过滤/合并**：`bcftools`
- **可视化**：`python3` + `matplotlib`
- **脚本运行**：`bash`，建议 `conda` 环境 `snpcaller`

```bash
conda create -n snpcaller -c bioconda bwa samtools bedtools bcftools gatk
conda activate snpcaller
pip install matplotlib
```

## 结果说明（含"表格错位"修复说明）

- **位点集**：过滤后 SNP 数量依输入样本与区域而定，由 pipeline 自动统计并打印（见 `gatk_call_summary.tsv`）。
- **基因型矩阵**：`genotypes_gatk.txt` 为合并多样本后的 `%CHROM %POS %REF %ALT %GT` 表。
- **位点对照表**：`results/gatk/upstream/lollipop/variant_positions.csv` 为**逗号分隔、UTF-8 BOM** 编码的 SNP 对照表
  （位置/变化/DP/AF/QUAL/各样本 GT/携带样本数）。
  - **为什么不会再"错位"**：早期版本导出的是 `.tsv`（tab 分隔），中文 Excel 默认不会按 tab 拆列，导致整行粘在一个单元格里。
    本仓库版本**固定输出 `.csv`（逗号分隔）+ BOM**，与输入数据无关——无论代入什么样本/数据，跑出来都是规整分列的，Excel 双击即用。
- **下游派生**：双等位金标准与下采样召回率分析可基于此位点集（合并 VCF / genotypes）进一步处理。

---

*本仓库为 GATK 版本上游 pipeline（姊妹仓库 ABCB5-SNP 使用 bcftools 检测）。*
