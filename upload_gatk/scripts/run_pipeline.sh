#!/bin/bash
set -euo pipefail

# ============================================================
# ABCB5-SNP-GATK · 上游分析主流程
# ------------------------------------------------------------
# 比对 (BWA-MEM) → GATK HaplotypeCaller → 双等位 SNP 过滤
#   → 合并多样本 VCF → 基因型矩阵 → 覆盖度统计 → 坐标转换
# 样本与区域由 scripts/config.sh 决定，不写死任何具体样本。
# 下游分析（饱和召回、金标准提纯、可视化等）不在此仓库中。
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ---- 激活环境（bwa / samtools / bcftools / bedtools）----
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snpcaller

echo "=========================================="
echo "  ABCB5 SNP Analysis — GATK Upstream"
echo "  PROJECT_DIR = $PROJECT_DIR"
echo "=========================================="

# ---- 检查必要文件 ----
for f in "$REF_GENOME" "$BED_FILE" "$SAMPLE_INFO"; do
    [ -f "$f" ] || { echo "错误: 文件不存在 - $f"; exit 1; }
done
[ -d "$FASTQ_DIR" ] && [ -n "$(ls -A "$FASTQ_DIR" 2>/dev/null)" ] || \
    { echo "错误: FASTQ 目录为空或不存在 - $FASTQ_DIR"; exit 1; }

mkdir -p "$ALIGN_DIR" "$VARIANT_DIR" "$SUMMARY_DIR" "$QC_DIR" "$COVERAGE_DIR" "$TEMP_DIR"

# ============================================================
# 准备参考索引（bwa / samtools / gatk）
# ============================================================
echo "[init] 准备参考索引 ..."
[ -f "${REF_GENOME}.bwt" ] || bwa index "$REF_GENOME"
[ -f "${REF_GENOME}.fai" ] || samtools faidx "$REF_GENOME"
[ -f "${REF_GENOME}.dict" ] || "$GATK_BIN" CreateSequenceDictionary -R "$REF_GENOME"

# ---- 读取样本列表（从 CSV，跳过表头）----
SAMPLES=$(tail -n +2 "$SAMPLE_INFO" | cut -d',' -f1 | sed 's/\r//g' | tr '\n' ' ')
echo "样本列表: $SAMPLES"

# ============================================================
# Step 1: 比对
# ============================================================
echo "Step 1: 序列比对 (BWA-MEM + bedtools 限区域) ..."
for sample in $SAMPLES; do
    echo "  处理 $sample ..."
    R1="${FASTQ_DIR}/${sample}_R1.fq.gz"
    R2="${FASTQ_DIR}/${sample}_R2.fq.gz"
    [ -f "$R1" ] && [ -f "$R2" ] || { echo "    警告: 未找到 $sample 的 FASTQ，跳过"; continue; }
    bwa mem -t "$NTHREADS" -R "@RG\tID:$sample\tSM:$sample\tLB:lib1\tPL:ILLUMINA" \
        "$REF_GENOME" "$R1" "$R2" | \
        samtools view -bS - | \
        bedtools intersect -abam stdin -b "$BED_FILE" | \
        samtools sort -o "${ALIGN_DIR}/${sample}.sorted.bam"
    samtools index "${ALIGN_DIR}/${sample}.sorted.bam"
    samtools flagstat "${ALIGN_DIR}/${sample}.sorted.bam" > "${QC_DIR}/${sample}.flagstat"
    echo "    完成 $sample"
done

# ============================================================
# Step 2: GATK HaplotypeCaller 变异检测 + 双等位 SNP 过滤
# ============================================================
echo "Step 2: GATK HaplotypeCaller ..."
call_gatk() {  # $1=inbam $2=outrawvcf
    "$GATK_BIN" --java-options "-Xmx${GATK_MEM}" HaplotypeCaller \
        -R "$REF_GENOME" -I "$1" -O "$2" \
        --native-pair-hmm-threads "$NTHREADS" --sample-ploidy "$SAMPLE_PLOIDY" \
        >/dev/null 2>>"$TEMP_DIR/gatk_errors.log"
}

# ---- 双等位 SNP 过滤口径（与已发布 p3 结果一致）----
FILTER="MIN(FMT/DP)>=$MIN_DP && MQ>=$MIN_MQ && (FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1]))>=$MIN_VAF && (GT=\"0/1\" || GT=\"1/1\")"

SUMMARY="$VARIANT_DIR/gatk_call_summary.tsv"
echo -e "SampleID\tRawSNPs\tFilteredSNPs" > "$SUMMARY"
VCFS=""
for sample in $SAMPLES; do
    BAM="${ALIGN_DIR}/${sample}.sorted.bam"
    [ -f "$BAM" ] || { echo "  跳过 $sample (无 BAM)"; continue; }
    RAW="$VARIANT_DIR/${sample}.gatk.raw.vcf.gz"
    FILT="$VARIANT_DIR/${sample}.gatk.f.vcf.gz"
    call_gatk "$BAM" "$RAW"
    bcftools view -v snps -i "$FILTER" "$RAW" -Oz -o "$FILT"
    bcftools index "$FILT" 2>>"$TEMP_DIR/gatk_errors.log"
    RAW_N=$(bcftools view -H "$RAW" | wc -l)
    FILT_N=$(bcftools view -H "$FILT" | wc -l)
    echo "=== $sample : GATK raw=$RAW_N 过滤后=$FILT_N ==="
    echo -e "${sample}\t${RAW_N}\t${FILT_N}" >> "$SUMMARY"
    VCFS="${VCFS:+$VCFS }$FILT"
    rm -f "$RAW"
done

# ============================================================
# Step 3: 合并多样本 VCF + 基因型矩阵
# ============================================================
echo "Step 3: 合并 VCF 并生成基因型矩阵 ..."
COMB="$VARIANT_DIR/filtered_snps.combined.vcf.gz"
GENO="$SUMMARY_DIR/genotypes.txt"
if [ "$(echo $VCFS | wc -w)" -gt 1 ]; then
    bcftools merge -Oz -o "$COMB" $VCFS 2>>"$TEMP_DIR/gatk_errors.log"
else
    cp $VCFS "$COMB" 2>>"$TEMP_DIR/gatk_errors.log"
fi
bcftools index "$COMB" 2>>"$TEMP_DIR/gatk_errors.log"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$COMB" > "$GENO"
HEADER="CHROM\tPOS\tREF\tALT"
for sample in $SAMPLES; do HEADER="$HEADER\t$sample"; done
sed -i "1i $HEADER" "$GENO"
echo "  基因型矩阵: $GENO"

# ============================================================
# Step 4: 覆盖度统计
# ============================================================
echo "Step 4: 覆盖度统计 ..."
TOTAL_BASES=$(awk '{sum+=$3-$2} END {print sum}' "$BED_FILE")
echo -e "Sample\tMean_depth\tCovered_bases\tTotal_bases\tCoverage_pct" > "${COVERAGE_DIR}/summary.tsv"
for sample in $SAMPLES; do
    depth_info=$(samtools depth -b "$BED_FILE" "${ALIGN_DIR}/${sample}.sorted.bam" 2>/dev/null | \
        awk '{sum+=$3; count++} END {if(count>0) printf "%s\t%d\t%.2f", sum/count, count, (count*100)/ENVIRON["TOTAL"]; else print "0\t0\t0.00"}' TOTAL="$TOTAL_BASES")
    echo -e "$sample\t$depth_info" >> "${COVERAGE_DIR}/summary.tsv"
done
cat "${COVERAGE_DIR}/summary.tsv"

# ============================================================
# Step 5: 坐标转换（相对 → 绝对）
# ============================================================
echo "Step 5: 坐标转换（相对 → 绝对）..."
INPUT_FILE="$GENO"
OUTPUT_FILE="${SUMMARY_DIR}/genotypes_abs.txt"
head -1 "$INPUT_FILE" > "$OUTPUT_FILE"
tail -n +2 "$INPUT_FILE" | while read line; do
    rel_pos=$(echo "$line" | cut -f2)
    if [[ "$rel_pos" =~ ^[0-9]+$ ]]; then
        abs_pos=$((rel_pos + OFFSET - 1))
        echo "$line" | awk -v abs="$abs_pos" 'BEGIN{OFS="\t"} {$2=abs; print}'
    else
        echo "$line"
    fi
done >> "$OUTPUT_FILE"
echo "  绝对坐标: $OUTPUT_FILE"

echo "=========================================="
echo "  GATK Upstream Pipeline 完成！"
echo "=========================================="
echo "输出文件:"
echo "  - $VARIANT_DIR/filtered_snps.combined.vcf.gz"
echo "  - $SUMMARY_DIR/genotypes.txt"
echo "  - $SUMMARY_DIR/genotypes_abs.txt (绝对坐标，用于 NCBI 查询)"
echo "  - $COVERAGE_DIR/summary.tsv"
echo "  - $QC_DIR/*.flagstat"
echo "  - $TEMP_DIR/gatk_errors.log"
echo "=========================================="
