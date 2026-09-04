#!/bin/bash
# ============================================================
# ABCB5-SNP-GATK · 上游 pipeline 配置
# ------------------------------------------------------------
# 所有路径与阈值均外部化，便于将本仓库嵌入任意数据目录。
# 修改本文件即可适配你的参考基因组、目标区域、样本与过滤阈值。
# ============================================================

# 自动定位仓库根目录（本文件位于 <repo>/scripts/config.sh）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# ---- 参考与区域（用户自备）----
REF_GENOME="${PROJECT_DIR}/data/reference.fa"        # 参考基因组（如 ABCB5 截取区域）
BED_FILE="${PROJECT_DIR}/data/targets.bed"            # 目标区域 BED
SAMPLE_INFO="${PROJECT_DIR}/config/sample_info.csv"   # 样本清单：sample_name,group
FASTQ_DIR="${PROJECT_DIR}/input"                      # <样本名>_R1.fq.gz / _R2.fq.gz

# ---- GATK 可执行文件 ----
# 默认使用本地 GATK 安装；可在环境变量中覆盖，例如：
#   GATK_BIN=/path/to/gatk-4.6.0.0/gatk bash scripts/run_pipeline.sh
GATK_BIN="${GATK_BIN:-/home/huangzy/gatk/gatk-4.6.0.0/gatk}"

# ---- 过滤阈值（双等位 SNP 口径，与已发布 p3 结果一致）----
MIN_DP=30        # 每个样本最小总深度
MIN_MQ=30        # 位点最小 Mapping Quality
MIN_VAF=0.30     # 杂合/纯合均要求 ALT 等位基因 VAF 不低于此值

# ---- GATK 运行参数 ----
SAMPLE_PLOIDY=2
GATK_MEM="8g"
NTHREADS=4

# ---- 输出目录（均为运行产物，已列入 .gitignore）----
OUTPUT_DIR="${PROJECT_DIR}/results"
ALIGN_DIR="${OUTPUT_DIR}/alignment"
VARIANT_DIR="${OUTPUT_DIR}/variants"
SUMMARY_DIR="${OUTPUT_DIR}/summary"
QC_DIR="${OUTPUT_DIR}/qc"
COVERAGE_DIR="${OUTPUT_DIR}/coverage"
TEMP_DIR="${OUTPUT_DIR}/temp"

# ---- 绝对坐标偏移（将相对坐标转为基因组绝对位置，依你的参考修改）----
OFFSET=20594465
