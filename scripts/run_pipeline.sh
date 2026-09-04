#!/bin/bash
# ============================================================
# ABCB5-SNP-GATK 上游 pipeline (GATK HaplotypeCaller)
#   建索引 → (比对) → HaplotypeCaller → 双等位过滤 → 合并 → 基因型矩阵 → (lollipop 图)
# 所有配置见同目录 config.sh；不写死任何样本或结果。
# 用法:
#   bash scripts/run_pipeline.sh
#   SKIP_ALIGN=true bash scripts/run_pipeline.sh        # 已有 BAM 时跳过比对
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

# 激活 conda（可选）
if [ -n "${CONDA_ENV:-}" ]; then
  source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
  conda activate "$CONDA_ENV" 2>/dev/null || true
fi

echo "=========================================="
echo "  ABCB5 SNP Analysis Pipeline (GATK)"
echo "  PROJECT_DIR = $PROJECT_DIR"
echo "=========================================="

# 检查必要文件
for f in "$REF_GENOME" "$BED_FILE" "$SAMPLE_INFO"; do
  [ -f "$f" ] || { echo "错误: 文件不存在 - $f"; exit 1; }
done
if [ "$SKIP_ALIGN" != "true" ]; then
  if [ ! -d "$FASTQ_DIR" ] || [ -z "$(ls -A "$FASTQ_DIR" 2>/dev/null)" ]; then
    echo "错误: FASTQ 目录为空或不存在 - $FASTQ_DIR (若已有 BAM, 请设 SKIP_ALIGN=true)"
    exit 1
  fi
fi

mkdir -p "$ALIGN_DIR" "$VCF_DIR" "$LOLLIPOP_DIR"

# 读取样本列表（从 CSV 读取，跳过表头与注释行）
SAMPLES=$(tail -n +2 "$SAMPLE_INFO" | cut -d',' -f1 | sed 's/\r//g' | grep -v '^#' | grep -v '^[[:space:]]*$' | tr '\n' ' ')
echo "样本列表: $SAMPLES"

# ---------- 建索引 ----------
echo "== 构建参考索引 =="
[ -f "${REF_GENOME}.bwt" ] || bwa index "$REF_GENOME"
[ -f "${REF_GENOME}.fai" ] || samtools faidx "$REF_GENOME"
DICT="${REF_GENOME%.*}.dict"
[ -f "$DICT" ] || "$GATK" CreateSequenceDictionary -R "$REF_GENOME"

# ---------- 比对（可选）----------
if [ "$SKIP_ALIGN" != "true" ]; then
  echo "== Step 1: 比对 (BWA-MEM) =="
  for sample in $SAMPLES; do
    BAM="${ALIGN_DIR}/${sample}.sorted.bam"
    [ -f "$BAM" ] && { echo "  $sample 已有 BAM, 跳过"; continue; }
    R1="${FASTQ_DIR}/${sample}_R1.fq.gz"; R2="${FASTQ_DIR}/${sample}_R2.fq.gz"
    [ -f "$R1" ] && [ -f "$R2" ] || { echo "  警告: 未找到 $sample FASTQ, 跳过"; continue; }
    bwa mem -t 4 -R "@RG\tID:$sample\tSM:$sample\tLB:lib1\tPL:ILLUMINA" \
      "$REF_GENOME" "$R1" "$R2" | samtools view -bS - | \
      bedtools intersect -abam stdin -b "$BED_FILE" | samtools sort -o "$BAM"
    samtools index "$BAM"
    echo "  完成 $sample"
  done
fi

# ---------- HaplotypeCaller + 双等位过滤 ----------
echo "== Step 2: HaplotypeCaller + 双等位过滤 =="
SUMMARY="${VCF_DIR}/gatk_call_summary.tsv"
echo -e "SampleID\tRawSNPs\tFilteredSNPs" > "$SUMMARY"
VCFS=""
for sample in $SAMPLES; do
  BAM="${ALIGN_DIR}/${sample}.sorted.bam"
  [ -f "$BAM" ] || { echo "skip $sample (no bam)"; continue; }
  RAW="${VCF_DIR}/${sample}.gatk.raw.vcf.gz"
  FILT="${VCF_DIR}/${sample}.gatk.f.vcf.gz"
  "$GATK" --java-options "-Xmx8g" HaplotypeCaller -R "$REF_GENOME" -I "$BAM" -O "$RAW" \
    --native-pair-hmm-threads 4 --sample-ploidy "$PLOIDY" 2>>"${VCF_DIR}/errors.log"
  bcftools view -v snps -i "$FILTER_EXPR" "$RAW" -Oz -o "$FILT"
  bcftools index "$FILT" 2>>"${VCF_DIR}/errors.log"
  RAW_N=$(bcftools view -H "$RAW" | wc -l)
  FILT_N=$(bcftools view -H "$FILT" | wc -l)
  echo "=== $sample : GATK raw=$RAW_N 过滤后=$FILT_N ==="
  echo -e "${sample}\t${RAW_N}\t${FILT_N}" >> "$SUMMARY"
  VCFS="${VCFS:+$VCFS }$FILT"
  rm -f "$RAW"
done

# ---------- 合并 + 基因型矩阵 ----------
echo "== Step 3: 合并 + 基因型矩阵 =="
COMB="${VCF_DIR}/gatk.combined.f.vcf.gz"
GENO="${VCF_DIR}/genotypes_gatk.txt"
N_VCF=$(echo $VCFS | wc -w)
if [ "$N_VCF" -gt 1 ]; then
  bcftools merge -Oz -o "$COMB" $VCFS 2>>"${VCF_DIR}/errors.log"
elif [ "$N_VCF" -eq 1 ]; then
  cp $VCFS "$COMB" 2>>"${VCF_DIR}/errors.log"
else
  echo "错误: 没有任何过滤后 VCF 可合并"; exit 1
fi
bcftools index "$COMB" 2>>"${VCF_DIR}/errors.log"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$COMB" > "$GENO"
echo "合并 VCF -> $COMB ($(wc -l < "$GENO") 个位点)"

# ---------- lollipop 图（可选）----------
if [ -f "$DOMAINS_FILE" ]; then
  echo "== Step 4: lollipop 图 + 位点对照表 (variant_positions.csv) =="
  export GATK_COMBINED_VCF="$COMB"
  export DOMAINS_FILE_PATH="$DOMAINS_FILE"
  export LOLLIPOP_DIR="$LOLLIPOP_DIR"
  python3 "$SCRIPT_DIR/plot_lollipop_gatk.py"
else
  echo "== 跳过 lollipop 图（缺少 config/domains.txt）=="
fi

echo "=========================================="
echo "  Pipeline 完成！"
echo "  产物目录: $VCF_DIR"
echo "    - gatk.combined.f.vcf.gz   合并多样本 VCF"
echo "    - genotypes_gatk.txt       基因型矩阵"
echo "    - gatk_call_summary.tsv    各样本原始/过滤 SNP 数"
echo "    - lollipop/                各样本 SNP 分布图 + variant_positions.csv"
echo "=========================================="
