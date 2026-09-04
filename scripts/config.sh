#!/bin/bash
# ============================================================
# ABCB5-SNP-GATK 配置 —— 所有路径/阈值外部化，不写死任何样本或结果。
# 本文件被 scripts/run_pipeline.sh 读取；
# run_pipeline.sh 会把关键路径以环境变量形式传给 plot_lollipop_gatk.py。
# ============================================================

# 仓库根目录：自动定位本脚本所在位置的上一级，便于嵌入任意数据目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# ---- 参考与区域 ----
REF_GENOME="${PROJECT_DIR}/data/reference.fa"
BED_FILE="${PROJECT_DIR}/data/targets.bed"

# ---- 输入样本清单 (CSV: sample_name,group) ----
SAMPLE_INFO="${PROJECT_DIR}/config/sample_info.csv"

# ---- 输入测序数据 (FASTQ) ----
FASTQ_DIR="${PROJECT_DIR}/input"

# ---- GATK 可执行文件（留空则依赖 PATH 中的 gatk）----
GATK="${GATK:-/home/huangzy/gatk/gatk-4.6.0.0/gatk}"

# ---- 变异过滤阈值（双等位口径，与 bcftools 姊妹仓库口径一致）----
MIN_DP=30
MIN_MQ=30
MIN_VAF=0.30
PLOIDY=2
# 过滤表达式：双等位 SNP，DP>=MIN_DP, MQ>=MIN_MQ, VAF>=MIN_VAF, GT 为 0/1 或 1/1
FILTER_EXPR='MIN(FMT/DP)>='"$MIN_DP"' && MQ>='"$MIN_MQ"' && (FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1]))>='"$MIN_VAF"' && (GT="0/1" || GT="1/1")'

# ---- 输出目录 ----
OUTPUT_DIR="${PROJECT_DIR}/results"
ALIGN_DIR="${OUTPUT_DIR}/alignment"
VCF_DIR="${OUTPUT_DIR}/gatk/upstream"        # 各样本 raw/filtered vcf + combined + genotypes
LOLLIPOP_DIR="${OUTPUT_DIR}/gatk/upstream/lollipop"

# ---- 可选：已有比对好的 BAM，则跳过比对步骤 ----
#   设为 true 时，pipeline 直接读取 ${ALIGN_DIR}/<样本>.sorted.bam
SKIP_ALIGN="${SKIP_ALIGN:-false}"

# ---- 可选：lollipop 图所需蛋白结构域文件（缺省则跳过绘图步骤）----
DOMAINS_FILE="${PROJECT_DIR}/config/domains.txt"

# ---- conda 环境（留空则不激活）----
CONDA_ENV="${CONDA_ENV:-snpcaller}"
