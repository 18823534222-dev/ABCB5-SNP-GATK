#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基于 GATK 上游 combined VCF 生成各样本 SNP 位点 lollipop 图（真实基因组坐标）+ 位点对照表。
- x 轴 = 真实基因组 POS
- 同一突变簇内的标签按垂直方向错层展开（扇形引线），避免重叠
- 结构域按真实坐标换算（domains.txt 坐标为 POS/100，乘以 100）
- 输出 variant_positions.csv：每个 SNP 的真实位置 / 变化 / 深度 / 频率 / 各样本 GT
  （逗号分隔 + UTF-8 BOM，Excel 双击即可正确分列，避免 tab 版在中文 Excel 里粘成一格）

依赖环境变量（由 scripts/run_pipeline.sh 注入）：
  GATK_COMBINED_VCF   合并后的多样本过滤 VCF
  DOMAINS_FILE_PATH   结构域坐标文件 (name<TAB>start<TAB>end, 坐标为 POS/100)
  LOLLIPOP_DIR        输出目录
"""
import os
import sys
import gzip
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.ticker as ticker

COMBINED = os.environ.get("GATK_COMBINED_VCF")
DOMAINS  = os.environ.get("DOMAINS_FILE_PATH")
OUTDIR   = os.environ.get("LOLLIPOP_DIR", "lollipop")

for name, val in (("GATK_COMBINED_VCF", COMBINED), ("DOMAINS_FILE_PATH", DOMAINS)):
    if not val:
        sys.exit(f"缺少环境变量 {name}（请通过 run_pipeline.sh 调用，或手动 export 后运行）")
os.makedirs(OUTDIR, exist_ok=True)

HET_COLOR = "#0072B2"   # 杂合
HOM_COLOR = "#D55E00"   # 纯合

plt.rcParams.update({
    "font.sans-serif": ["Microsoft YaHei", "Noto Sans SC", "SimSun", "DejaVu Sans"],
    "axes.unicode_minus": False,
})

# ---------- 读结构域 (坐标为 POS/100) ----------
domains = []
with open(DOMAINS) as f:
    next(f, None)
    for line in f:
        line = line.strip()
        if not line:
            continue
        p = line.split("\t")
        if len(p) < 3:
            continue
        domains.append((p[0], float(p[1]), float(p[2])))

def gt_type(gt):
    g = gt.replace("|", "/")
    if g in ("", ".", "./.", "0/0"):
        return None
    return "het" if ("/0" in g or "0/" in g) else "hom"

def parse_info(info):
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d

# ---------- 解析 combined VCF ----------
records = []
samples = []
sample_vars = { }
with gzip.open(COMBINED, "rt") as f:
    for line in f:
        if line.startswith("#CHROM"):
            cols = line.rstrip("\n").split("\t")
            samples = cols[9:]
            for s in samples:
                sample_vars[s] = []
            continue
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        chrom, pos = fields[0], int(fields[1])
        ref, alt, qual = fields[3], fields[4], fields[5]
        info = parse_info(fields[7])
        fmt = fields[8].split(":")
        gt_idx = fmt.index("GT") if "GT" in fmt else 0
        dp, af = info.get("DP", ""), info.get("AF", "")
        rec = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
               "qual": qual, "dp": dp, "af": af, "gt": {}}
        for i, s in enumerate(samples):
            parts = fields[9 + i].split(":")
            gt = parts[gt_idx] if gt_idx < len(parts) else fields[9 + i]
            rec["gt"][s] = gt
            t = gt_type(gt)
            if t is not None:
                sample_vars[s].append((pos, ref, alt, t))
        records.append(rec)

# ---------- 输出位点对照表（逗号分隔 CSV，Excel 直接分列）----------
#  注意：固定用逗号分隔 + UTF-8 BOM，与输入数据无关 —— 任何数据跑出来都是规整分列的。
csv_path = os.path.join(OUTDIR, "variant_positions.csv")
header = ["位置(POS)", "参考碱基", "变异碱基", "变化", "深度DP", "等位频率AF",
          "QUAL"] + samples + ["携带样本数"]
with open(csv_path, "w", newline="", encoding="utf-8-sig") as fo:
    w = csv.writer(fo, delimiter=",")
    w.writerow(header)
    for rec in sorted(records, key=lambda r: r["pos"]):
        n_car = sum(1 for s in samples if gt_type(rec["gt"][s]) is not None)
        row = [rec["pos"], rec["ref"], rec["alt"], f'{rec["ref"]}>{rec["alt"]}',
               rec["dp"], rec["af"], rec["qual"]]
        row += [rec["gt"][s] for s in samples]
        row += [n_car]
        w.writerow(row)
print("位点对照表已写出:", csv_path, f"({len(records)} 个变异)")

# ---------- 画每个样本（真实坐标 + 簇内标签错层展开）----------
all_pos = [r["pos"] for r in records]
min_pos, max_pos = min(all_pos), max(all_pos)
x_lo, x_hi = max(0, min_pos - 5000), max_pos + 8000
domain_max = max((d[2] for d in domains), default=0) * 100
protein_hi = max(max_pos, domain_max) + 8000

CLUSTER_GAP = 10000.0  # 相邻变异间距 > 此值则另起一簇（保证簇间标签不横撞）
Y0 = 0.98              # 第一层标签 y
YSTEP = 0.15           # 标签层间距

for s in samples:
    vars_ = sorted(sample_vars[s])
    # 分簇
    clusters, cur = [], []
    for v in vars_:
        if cur and (v[0] - cur[-1][0]) > CLUSTER_GAP:
            clusters.append(cur)
            cur = []
        cur.append(v)
    if cur:
        clusters.append(cur)
    max_levels = max((len(c) for c in clusters), default=1)
    y_top = Y0 + (max_levels - 1) * YSTEP + 0.35

    fig, ax = plt.subplots(figsize=(11, 3.8 + 0.18 * max_levels))
    ax.add_patch(Rectangle((0, -0.12), protein_hi, 0.24,
                           facecolor="#E8E8E8", edgecolor="#999999", lw=0.8))
    palette = ["#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462"]
    for i, (name, st, en) in enumerate(domains):
        a, b = st * 100, en * 100
        ax.add_patch(Rectangle((a, -0.12), b - a, 0.24,
                               facecolor=palette[i % len(palette)], alpha=0.7, edgecolor="none"))
        ax.text((a + b) / 2, 0.02, name, ha="center", va="bottom", fontsize=7, color="#333333")
    top = 0.85
    for cl in clusters:
        cx = sum(p for p, _, _, _ in cl) / len(cl)   # 簇中心 x（标签堆叠处）
        for i, (pos, ref, alt, t) in enumerate(cl):
            c = HET_COLOR if t == "het" else HOM_COLOR
            ax.plot([pos, pos], [0, top], color=c, lw=1.0, zorder=3)
            ax.scatter([pos], [top], s=28, color=c, edgecolor="white", lw=0.6, zorder=4)
            ly = Y0 + i * YSTEP
            ax.annotate(f"{pos:,} {ref}>{alt}", xy=(pos, top), xytext=(cx, ly),
                        ha="center", va="bottom", fontsize=6.0, color=c,
                        bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.75),
                        arrowprops=dict(arrowstyle="-", color=c, lw=0.4, ls="--"))
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(-0.3, y_top)
    ax.set_yticks([])
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.set_xticks(range(0, int(protein_hi) + 1, 20000))
    ax.set_xlabel("Genomic Position (ABCB5, bp)")
    ax.set_title(f"{s} - SNP 位点分布（杂合蓝 / 纯合朱红）", fontsize=11, loc="left")
    ax.grid(axis="x", color="#EEEEEE", lw=0.6)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="both", which="both", direction="in")
    fig.tight_layout()
    png = os.path.join(OUTDIR, f"{s}_lollipop.png")
    pdf = os.path.join(OUTDIR, f"{s}_lollipop.pdf")
    fig.savefig(png, dpi=200)
    fig.savefig(pdf, dpi=200)
    plt.close(fig)
    print(f"{s}: {len(vars_)} 个变异位点 -> {os.path.basename(png)}")
