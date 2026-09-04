# ABCB5 SNP Analysis — GATK Upstream Pipeline

> 一套**可复现、可嵌入任意数据**的 germline SNP 上游检测与分型流程，基于 **GATK HaplotypeCaller**。
> 从原始测序数据（FASTQ）到高质量 SNP 位点表、基因型矩阵与覆盖度统计，
> 样本与参考区域均由配置文件决定，不锁定任何具体样本或结果。
>
> **本仓库仅包含上游 pipeline**（比对 → 变异检测 → 过滤 → 基因型 → 覆盖度）。
> 下游分析（金标准提纯、下采样召回率、可视化等）不在此仓库中。
> 变异检测工具为 **GATK**；bcftools 版本见姊妹仓库 **ABCB5-SNP-bcftools**。

---

## 目录

- [概览](#概览)
- [分析流程](#分析流程)
- [样本](#样本)
- [文件结构](#文件结构)
- [快速使用](#快速使用)
- [依赖环境](#依赖环境)
- [结果说明](#结果说明)

---

## 概览

本仓库交付一套**完整、可复现**的上游分析 pipeline，用于从原始测序数据
（FASTQ）中检出 ABCB5 区域内的 germline SNP，并产出：

- 过滤后的高质量 SNP 位点集（合并 VCF + 易读表格）
- 各样本基因型矩阵（相对坐标 / 绝对坐标）
- 覆盖度与质量汇总统计

下游分析（如双等位金标准构建、下采样召回率分析等）可基于此位点集进一步派生。

## 分析流程

| 步骤 | 工具 | 说明 |
|------|------|------|
| 比对 | BWA-MEM | 基于用户提供的参考基因组与目标区域 BED，限制到目标区域 |
| 变异检测 | `GATK HaplotypeCaller` | 逐样本 call，输出原始 VCF |
| 过滤 | 双等位 SNP 口径 | `TYPE=snp`，`DP>=30`，`MQ>=30`，`VAF>=0.30`，`GT` 为 `0/1` 或 `1/1` |
| 合并 | `bcftools merge` | 合并为多样本 VCF |
| 基因型 | `bcftools query` | 输出基因型矩阵 |

> 过滤阈值为外部化参数，定义在 `scripts/config.sh` 中，可按需调整。

## 样本

样本由 `config/sample_info.csv` 决定，**不写死在代码中**。
在该 CSV 中按 `sample_name,group` 格式逐行列出你的样本即可；
样本名须与 `input/` 目录下 `<样本名>_R1.fq.gz` / `<样本名>_R2.fq.gz` 对应。
（仓库自带一份示例样本清单，仅作模板，可替换为任意数据。）

## 文件结构

```
ABCB5-SNP-GATK/
├── README.md
├── .gitignore
├── scripts/
│   ├── run_pipeline.sh          # 上游 pipeline 主脚本（一键复现）
│   └── config.sh                # 所有路径/阈值配置（外部化）
├── config/
│   └── sample_info.csv          # 样本清单（用户自定义）
├── data/                        # 用户自备：reference.fa / targets.bed
├── input/                       # 用户自备：<样本名>_R1.fq.gz / _R2.fq.gz
└── results/                     # 运行产物（不纳入版本控制）
    ├── alignment/               # 比对 BAM
    ├── variants/                # 各样本 filtered VCF + 合并 VCF + 统计
    ├── summary/                 # 基因型矩阵（相对/绝对坐标）
    ├── qc/                      # flagstat 比对率
    └── coverage/                # 覆盖度统计
```

> 运行产物（`results/` 下各目录）由 pipeline 动态生成，已列入 `.gitignore`，
> 不会随仓库分发。

## 快速使用

```bash
# 1. 克隆仓库
git clone https://github.com/18823534222-dev/ABCB5-SNP-GATK.git
cd ABCB5-SNP-GATK

# 2. 激活 conda 环境（含 bwa / samtools / bcftools / bedtools）
conda activate snpcaller

# 3. 配置 GATK 可执行文件路径（如非默认路径，用环境变量覆盖）
export GATK_BIN=/path/to/gatk-4.6.0.0/gatk

# 4. 运行上游 pipeline（需自备输入 FASTQ 与参考）
bash scripts/run_pipeline.sh
```

> 注：仓库本身**不包含运行产物**（已列入 `.gitignore`）。
> 克隆后请自备 `data/`（参考与 BED）与 `input/`（样本 FASTQ），
> 运行 `bash scripts/run_pipeline.sh` 即可在本机生成 `results/` 下各结果。

## 依赖环境

- **比对**：`bwa` (BWA-MEM)
- **处理**：`samtools`
- **变异检测**：`GATK` (HaplotypeCaller, 4.6.0.0)
- **过滤/合并/基因型**：`bcftools`
- **区域限制**：`bedtools`
- **脚本运行**：`bash`，建议 `conda` 环境 `snpcaller`

```bash
conda create -n snpcaller -c bioconda bwa samtools bcftools bedtools
conda activate snpcaller
# GATK 需另行安装（见 https://github.com/broadinstitute/gatk），
# 并将其 gatk 可执行文件路径写入 scripts/config.sh 的 GATK_BIN。
```

## 结果说明

- **位点集**：过滤后 SNP 数量依输入样本与区域而定，由 pipeline 自动统计并打印。
- **绝对坐标**：`summary/genotypes_abs.txt` 使用脚本中 `OFFSET`（默认 `20594465`，
  对应示例目标区域）将相对坐标转换为绝对基因组位置，可直接用于 NCBI / UCSC 查询。
- **下游派生**：双等位金标准与下采样召回率分析可基于此位点集进一步处理。
- **与 bcftools 版本对比**：本仓库为 GATK 上游实现；方法学对比见姊妹仓库
  **ABCB5-SNP-bcftools** 及其下游对比脚本。

---

*本仓库为 GATK 上游版本（v1.0）。*
