#!/usr/bin/env python3
"""
单数据集流水线：
  1. fasterq-dump 解压 SRA -> fastq（数据集内样本多线程）
  2. fastp QC（数据集内样本多线程）
  3. HISAT2 比对（数据集内样本多线程）
  4. StringTie 定量（数据集内样本多线程）
  5. 删除 fastq 文件（释放磁盘空间）
  6. 写出完成标志文件

用法:
  python 02_dataset_pipeline.py \
      --dataset_id GSE119834 \
      --species homo \
      --sra_dir result/1_download_geo/homo/rawdata/GSE119834 \
      --fastq_dir result/1_download_geo/homo/rawdata/GSE119834 \
      --qc_dir result/2_QC/homo/GSE119834 \
      --align_dir result/3_Align_Filter/homo/GSE119834 \
      --threads 8 \
      --maxparallel 4 \
      --output_marker result/3_Align_Filter/homo/GSE119834/dataset_finished.txt
"""

import os
import sys
import argparse
import subprocess
import glob
import shutil
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime


def ts_log(msg: str):
    """带时间戳的日志输出，格式: [2026-08-08 10:00:01] msg"""
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

# ─────────────────────────────────────────────
# 项目内置注释路径（workflow/anno/）
# 所有路径相对于项目根目录（Snakefile 所在目录）
# 可通过 --anno_base 和 --hisat2_index 命令行参数覆盖
# ─────────────────────────────────────────────
_DEFAULT_ANNO_BASE = "workflow/anno"

# HISAT2 索引前缀（不含 .ht2 后缀）
# 注意：与实际目录名大小写保持一致（humanhisat2Index，全小写）
_DEFAULT_HISAT2_INDEX = {
    "homo":  f"{_DEFAULT_ANNO_BASE}/homo/humanhisat2Index/GRCh38",
    "mouse": f"{_DEFAULT_ANNO_BASE}/mouse/mouseHisat2Index/GRCm38",
}


def _build_annotations(anno_base):
    """根据 anno_base 构建 8 类注释的 GTF 文件与对应定量子目录"""
    return {
        "homo": [
            # mRNA（需用户放置 GENCODE v44 mRNA GTF）
            {"gff": f"{anno_base}/homo/mRNAanno/gencode.v44.mRNA.annotation.gtf",
             "dir": "mRNA/genecode/stringtie"},
            # eRNA
            {"gff": f"{anno_base}/homo/ncRNAanno/EnhancerAtlasv2.0_eRNA.hg38.gtf",
             "dir": "eRNA/EnhancerAtlas/stringtie"},
            {"gff": f"{anno_base}/homo/ncRNAanno/Ensemblv110_eRNA.gtf",
             "dir": "eRNA/Ensembl/stringtie"},
            {"gff": f"{anno_base}/homo/ncRNAanno/FANTOM5_eRNA.gtf",
             "dir": "eRNA/FANTOM5/stringtie"},
            # lncRNA
            {"gff": f"{anno_base}/homo/ncRNAanno/GENCODEv44_lncRNA.gtf",
             "dir": "lncRNA/GENCODE/stringtie"},
            {"gff": f"{anno_base}/homo/ncRNAanno/NONCODEv6_lncRNA.gtf",
             "dir": "lncRNA/NONCODE/stringtie"},
            # miRNA
            {"gff": f"{anno_base}/homo/ncRNAanno/miRBasev22.1_miRNA.gtf",
             "dir": "miRNA/miRBase/stringtie"},
            {"gff": f"{anno_base}/homo/ncRNAanno/MirGeneDBv2.1_miRNA.gtf",
             "dir": "miRNA/MirGeneDB/stringtie"},
        ],
        "mouse": [
            {"gff": f"{anno_base}/mouse/mRNAanno/gencode.vM25.mRNA.annotation.gtf",
             "dir": "mRNA/genecode/stringtie"},
            {"gff": f"{anno_base}/mouse/ncRNAanno/EnhancerAtlasv2.0_eRNA.gtf",
             "dir": "eRNA/EnhancerAtlas/stringtie"},
            {"gff": f"{anno_base}/mouse/ncRNAanno/Ensemblv110_eRNA.gtf",
             "dir": "eRNA/Ensembl/stringtie"},
            {"gff": f"{anno_base}/mouse/ncRNAanno/FANTOM5_eRNA.gtf",
             "dir": "eRNA/FANTOM5/stringtie"},
            {"gff": f"{anno_base}/mouse/ncRNAanno/GENCODEvM25_lncRNA.gtf",
             "dir": "lncRNA/GENCODE/stringtie"},
            {"gff": f"{anno_base}/mouse/ncRNAanno/NONCODEv6_lncRNA.gtf",
             "dir": "lncRNA/NONCODE/stringtie"},
            {"gff": f"{anno_base}/mouse/ncRNAanno/miRBasev22.1_miRNA.gtf",
             "dir": "miRNA/miRBase/stringtie"},
            {"gff": f"{anno_base}/mouse/ncRNAanno/MirGeneDBv2.1_miRNA.gtf",
             "dir": "miRNA/MirGeneDB/stringtie"},
        ],
    }


# 兼容旧代码：模块级变量（后续被 main() 中 args 覆盖）
HISAT2_INDEX = _DEFAULT_HISAT2_INDEX.copy()
ANNOTATIONS = _build_annotations(_DEFAULT_ANNO_BASE)

# PICARD_JAR 默认值：优先使用项目内置路径，可通过 --picard_jar 覆盖
_DEFAULT_PICARD = "workflow/env/picard-2.18.2/picard.jar"


def parse_args():
    parser = argparse.ArgumentParser(description="Per-dataset pipeline: decompress -> QC -> align -> quant -> cleanup")
    parser.add_argument("--dataset_id",    required=True, help="Dataset ID, e.g. GSE119834")
    parser.add_argument("--species",       required=True, help="Species: homo or mouse")
    parser.add_argument("--sra_dir",       required=True, help="Directory containing .sra files")
    parser.add_argument("--fastq_dir",     required=True, help="Directory to write raw fastq files (temp)")
    parser.add_argument("--qc_dir",        required=True, help="Output directory for clean fastq and QC reports")
    parser.add_argument("--align_dir",     required=True, help="Output directory for alignment and quantification")
    parser.add_argument("--threads",       type=int, default=8,  help="Threads per sample task")
    parser.add_argument("--maxparallel",   type=int, default=4,  help="Max parallel samples within dataset")
    parser.add_argument("--output_marker", required=True, help="Path to write finished marker file")
    parser.add_argument("--picard_jar",    default=_DEFAULT_PICARD,
                        help=f"Path to picard.jar (default: {_DEFAULT_PICARD})")
    parser.add_argument(
        "--strandedness",
        default="unstranded",
        choices=["unstranded", "forward", "reverse"],
        help=(
            "链特异性设置（对应 HISAT2 --rna-strandness 和 StringTie --rf/--fr）。\n"
            "  unstranded : 不设置链特异性参数（默认，兼容非链特异性文库）\n"
            "  forward    : FR / Illumina TruSeq Stranded mRNA (read1 正链)\n"
            "               → HISAT2: --rna-strandness FR；StringTie: --fr\n"
            "  reverse    : RF / dUTP 方法（最常见链特异性文库）\n"
            "               → HISAT2: --rna-strandness RF；StringTie: --rf\n"
            "建议通过 RSeQC infer_experiment.py 确认文库类型后再指定。"
        )
    )
    parser.add_argument(
        "--sample_chunk_size",
        type=int, default=0,
        help=(
            "滚动窗口大小：每次处理 x 个样本，处理完后立即删除 sra/fastq/bam，"
            "再处理下一批，防止磁盘占用膨胀。\n"
            "  0（默认）: 不分块，一次性处理所有样本（原始行为）\n"
            "  正整数   : 每批处理 x 个样本，推荐值 4-8"
        )
    )
    parser.add_argument(
        "--hisat2_index",
        default="",
        help=(
            "HISAT2 索引前缀路径（覆盖内置默认值）。\n"
            "例如: workflow/anno/homo/humanhisat2Index/GRCh38\n"
            "若不指定则使用内置默认值。"
        )
    )
    parser.add_argument(
        "--anno_base",
        default="",
        help=(
            "注释文件根目录（覆盖内置默认值 workflow/anno）。\n"
            "GTF 文件相对于此目录按 {species}/mRNAanno/ 等子目录组织。"
        )
    )
    parser.add_argument(
        "--remove_duplicates",
        default="true",
        choices=["true", "false"],
        help=(
            "Picard MarkDuplicates 是否删除重复读段。\n"
            "  true（默认）: REMOVE_DUPLICATES=true（删除 dup reads）\n"
            "  false       : REMOVE_DUPLICATES=false（仅标记，保留 dup reads）\n"
            "RNA-seq 定量通常建议 false（仅标记），但默认保持 true 以兼容现有结果。"
        )
    )
    return parser.parse_args()


# ─────────────────────────────────────────────
# Step 1: fasterq-dump
# ─────────────────────────────────────────────
def decompress_sra(sra_file, fastq_dir, threads):
    srr = os.path.splitext(os.path.basename(sra_file))[0]
    fq1 = os.path.join(fastq_dir, f"{srr}_1.fastq")
    fq2 = os.path.join(fastq_dir, f"{srr}_2.fastq")
    fq_se = os.path.join(fastq_dir, f"{srr}.fastq")

    if os.path.exists(fq1) or os.path.exists(fq_se):
        ts_log(f"[Skip-Decompress] {srr}: fastq 已存在")
        return True

    ts_log(f"[Decompress] {srr}: 开始解压 ...")
    t0 = datetime.now()
    # 使用 --split-3 替代 --split-files，更稳健地处理异常 paired 记录
    cmd = ["fasterq-dump", "--split-3", "--threads", str(threads), "-O", fastq_dir, sra_file]
    try:
        subprocess.run(cmd, check=True)
        elapsed = (datetime.now() - t0).total_seconds()
        ts_log(f"[Decompress-Done] {srr}: {elapsed:.1f}s")
        return True
    except subprocess.CalledProcessError as e:
        ts_log(f"[Decompress-Error] {srr}: {e}")
        return False


# ─────────────────────────────────────────────
# Step 2: fastp QC
# ─────────────────────────────────────────────
def run_fastp(srr, fastq_dir, qc_dir, threads):
    fq1 = os.path.join(fastq_dir, f"{srr}_1.fastq")
    fq2 = os.path.join(fastq_dir, f"{srr}_2.fastq")
    fq_se = os.path.join(fastq_dir, f"{srr}.fastq")

    out1 = os.path.join(qc_dir, f"{srr}_1.clean.fastq")
    out2 = os.path.join(qc_dir, f"{srr}_2.clean.fastq")
    out_se = os.path.join(qc_dir, f"{srr}.clean.fastq")

    report_dir = os.path.join(qc_dir, "qc_reports")
    os.makedirs(report_dir, exist_ok=True)
    html = os.path.join(report_dir, f"{srr}_fastp.html")
    json = os.path.join(report_dir, f"{srr}_fastp.json")

    # 跳过已完成
    if os.path.exists(out1) or os.path.exists(out_se):
        ts_log(f"[Skip-QC] {srr}: clean fastq 已存在")
        return True

    ts_log(f"[QC] {srr}: 开始 fastp ...")
    t0 = datetime.now()
    if os.path.exists(fq1) and os.path.exists(fq2):
        cmd = ["fastp", "-i", fq1, "-I", fq2, "-o", out1, "-O", out2,
               "-h", html, "-j", json, "-w", str(threads)]
    elif os.path.exists(fq1):
        cmd = ["fastp", "-i", fq1, "-o", out1, "-h", html, "-j", json, "-w", str(threads)]
    elif os.path.exists(fq_se):
        cmd = ["fastp", "-i", fq_se, "-o", out_se, "-h", html, "-j", json, "-w", str(threads)]
    else:
        ts_log(f"[QC-Error] {srr}: 找不到 fastq 文件")
        return False

    try:
        subprocess.run(cmd, check=True)
        elapsed = (datetime.now() - t0).total_seconds()
        ts_log(f"[QC-Done] {srr}: {elapsed:.1f}s")
        return True
    except subprocess.CalledProcessError as e:
        ts_log(f"[QC-Error] {srr}: {e}")
        return False


# ─────────────────────────────────────────────
# Step 3: HISAT2 比对
# ─────────────────────────────────────────────
def run_hisat2(srr, qc_dir, align_dir, species, threads, strandedness="unstranded", picard_jar=None, remove_duplicates=True):
    index = HISAT2_INDEX.get(species)
    if not index:
        ts_log(f"[Align-Error] {srr}: 未知物种 {species}")
        return False

    fq1 = os.path.join(qc_dir, f"{srr}_1.clean.fastq")
    fq2 = os.path.join(qc_dir, f"{srr}_2.clean.fastq")
    fq_se = os.path.join(qc_dir, f"{srr}.clean.fastq")

    hisat2_dir = os.path.join(align_dir, "hisat2file", srr)
    os.makedirs(hisat2_dir, exist_ok=True)

    dedup_bam = os.path.join(hisat2_dir, f"{srr}.dedup.bam")
    if os.path.exists(dedup_bam):
        ts_log(f"[Skip-Align] {srr}: dedup BAM 已存在")
        return True

    ts_log(f"[Align] {srr}: 开始 HISAT2 比对 (strandedness={strandedness}) ...")
    t0_align = datetime.now()
    sam_file  = os.path.join(hisat2_dir, "accepted_hits.sam")
    bam_file  = os.path.join(hisat2_dir, f"{srr}.sorted.bam")
    qc_log    = os.path.join(hisat2_dir, "QC_results.log")

    hisat2_cmd = [
        "hisat2", "-x", index, "-p", str(threads), "--dta",
        "--rg-id", srr, "--rg", f"SM:{srr}"
    ]

    # 链特异性：unstranded 不传参数，forward=FR，reverse=RF
    _strandmap = {"forward": "FR", "reverse": "RF"}
    if strandedness in _strandmap:
        hisat2_cmd += ["--rna-strandness", _strandmap[strandedness]]

    is_paired = os.path.exists(fq1) and os.path.exists(fq2)
    if is_paired:
        hisat2_cmd += ["-1", fq1, "-2", fq2]
    elif os.path.exists(fq1):
        hisat2_cmd += ["-U", fq1]
    elif os.path.exists(fq_se):
        hisat2_cmd += ["-U", fq_se]
    else:
        ts_log(f"[Align-Error] {srr}: 找不到 clean fastq")
        return False
    hisat2_cmd += ["-S", sam_file]

    try:
        with open(qc_log, "w") as log:
            subprocess.run(hisat2_cmd, check=True, stderr=log)

        subprocess.run(
            f"samtools view -bS {sam_file} | samtools sort -@ {threads} -o {bam_file}",
            shell=True, check=True)
        subprocess.run(f"samtools index {bam_file}", shell=True, check=True)

        if os.path.exists(sam_file):
            os.remove(sam_file)

        _picard = picard_jar or _DEFAULT_PICARD
        _remove_dup_flag = "true" if remove_duplicates else "false"
        dedup_cmd = [
            "java", "-Xmx15g", "-jar", _picard, "MarkDuplicates",
            f"I={bam_file}", f"O={dedup_bam}",
            f"METRICS_FILE={hisat2_dir}/{srr}.metrics",
            f"REMOVE_DUPLICATES={_remove_dup_flag}", "ASSUME_SORT_ORDER=coordinate"
        ]
        subprocess.run(dedup_cmd, check=True)
        subprocess.run(f"samtools index {dedup_bam}", shell=True, check=True)

        if os.path.exists(bam_file):
            os.remove(bam_file)
            bai = bam_file + ".bai"
            if os.path.exists(bai):
                os.remove(bai)

        elapsed = (datetime.now() - t0_align).total_seconds()
        ts_log(f"[Align-Done] {srr}: {elapsed:.1f}s")
        return True
    except subprocess.CalledProcessError as e:
        ts_log(f"[Align-Error] {srr}: {e}")
        return False


# ─────────────────────────────────────────────
# Step 4: StringTie 定量
# ─────────────────────────────────────────────
def run_stringtie(srr, align_dir, species, threads, strandedness="unstranded"):
    annotations = ANNOTATIONS.get(species, [])
    if not annotations:
        ts_log(f"[Quant-Error] {srr}: 未知物种 {species}")
        return False

    hisat2_dir = os.path.join(align_dir, "hisat2file", srr)
    dedup_bam  = os.path.join(hisat2_dir, f"{srr}.dedup.bam")

    if not os.path.exists(dedup_bam):
        ts_log(f"[Quant-Error] {srr}: dedup BAM 不存在，跳过定量")
        return False

    # StringTie 链特异性标志：unstranded 不加，forward→--fr，reverse→--rf
    _strand_flag = {"forward": ["--fr"], "reverse": ["--rf"]}.get(strandedness, [])

    ts_log(f"[Quant] {srr}: 开始 StringTie 定量 (strandedness={strandedness}) ...")
    t0_quant = datetime.now()
    for anno in annotations:
        out_dir   = os.path.join(align_dir, anno["dir"], srr)
        os.makedirs(out_dir, exist_ok=True)
        out_gtf   = os.path.join(out_dir, "transcripts.gtf")
        gene_abund = os.path.join(out_dir, "gene_abund.tab")

        if os.path.exists(out_gtf):
            continue

        cmd = (
            ["stringtie", "-p", str(threads), "-e", "-B"]
            + _strand_flag
            + ["-G", anno["gff"], "-A", gene_abund, "-o", out_gtf, dedup_bam]
        )
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            ts_log(f"[Quant-Error] {srr} ({anno['dir']}): {e}")

    elapsed = (datetime.now() - t0_quant).total_seconds()
    ts_log(f"[Quant-Done] {srr}: {elapsed:.1f}s")
    return True


# ─────────────────────────────────────────────
# Step 5: 删除 fastq 文件
# ─────────────────────────────────────────────
def cleanup_fastq(srr, fastq_dir, qc_dir):
    """删除原始 fastq 和 clean fastq，释放磁盘空间"""
    removed = []
    # 原始 fastq
    for pattern in [f"{srr}_1.fastq", f"{srr}_2.fastq", f"{srr}.fastq"]:
        p = os.path.join(fastq_dir, pattern)
        if os.path.exists(p):
            os.remove(p)
            removed.append(p)
    # clean fastq
    for pattern in [f"{srr}_1.clean.fastq", f"{srr}_2.clean.fastq", f"{srr}.clean.fastq"]:
        p = os.path.join(qc_dir, pattern)
        if os.path.exists(p):
            os.remove(p)
            removed.append(p)
    if removed:
        ts_log(f"[Cleanup] {srr}: 已删除 {len(removed)} 个 fastq 文件")
    return True


# ─────────────────────────────────────────────
# 单样本完整流水线
# ─────────────────────────────────────────────
def process_sample(srr, args):
    """对单个 SRR 执行完整流水线：解压 -> QC -> 比对 -> 定量 -> 清理"""
    print(f"\n{'='*50}")
    ts_log(f"[Pipeline] 开始处理样本: {srr} (数据集: {args.dataset_id})")
    print(f"{'='*50}")

    sra_file = os.path.join(args.sra_dir, f"{srr}.sra")
    if not os.path.exists(sra_file):
        ts_log(f"[Pipeline-Error] {srr}: SRA 文件不存在: {sra_file}")
        return False

    # 1. 解压
    if not decompress_sra(sra_file, args.fastq_dir, args.threads):
        ts_log(f"[Pipeline-Error] {srr}: 解压失败")
        return False

    # 2. QC
    if not run_fastp(srr, args.fastq_dir, args.qc_dir, args.threads):
        ts_log(f"[Pipeline-Error] {srr}: QC 失败")
        return False

    # 3. 比对（传入链特异性设置、picard 路径和去重标志）
    _remove_dup = getattr(args, "remove_duplicates", "true") == "true"
    if not run_hisat2(srr, args.qc_dir, args.align_dir, args.species,
                      args.threads,
                      getattr(args, "strandedness", "unstranded"),
                      getattr(args, "picard_jar", None),
                      remove_duplicates=_remove_dup):
        ts_log(f"[Pipeline-Error] {srr}: 比对失败")
        return False

    # 4. StringTie 定量（传入链特异性设置）
    if not run_stringtie(srr, args.align_dir, args.species,
                         args.threads, getattr(args, "strandedness", "unstranded")):
        ts_log(f"[Pipeline-Error] {srr}: 定量失败")
        return False

    # 5. 清理 fastq（释放空间）
    cleanup_fastq(srr, args.fastq_dir, args.qc_dir)

    ts_log(f"[Pipeline-Done] {srr} 全部完成")
    return True


# ─────────────────────────────────────────────
# 工具函数：滚动窗口清理
# ─────────────────────────────────────────────

def cleanup_chunk(success_srrs, args, keep_bam=False, keep_sra=False):
    """
    滚动窗口清理：只对**成功完成**的样本删除 fastq/bam，释放磁盘。
    失败的样本保留 .sra，以便 --rerun-incomplete 重试。
    只在 sample_chunk_size > 0 时调用（分块模式）。

    注意：.sra 文件不在此处删除——统一交给 rule shell 末尾的兜底清理
    （只有整个 GSE 成功时才执行）。这样强行停止/崩溃恢复后，
    已处理样本的 .sra 仍然存在，可以重试而非重新下载。
    """
    for srr in success_srrs:
        # .sra 保留：由 rule 末尾统一清理（防止进程崩溃后数据丢失）
        ts_log(f"[ChunkClean] {srr}: 保留 .sra（由 rule 末尾统一清理）")

        # 删除原始 fastq 和 clean fastq（process_sample 已调 cleanup_fastq，这里兜底）
        for pat in [f"{srr}_1.fastq", f"{srr}_2.fastq", f"{srr}.fastq"]:
            p = os.path.join(args.fastq_dir, pat)
            if os.path.exists(p):
                os.remove(p)
        for pat in [f"{srr}_1.clean.fastq", f"{srr}_2.clean.fastq", f"{srr}.clean.fastq"]:
            p = os.path.join(args.qc_dir, pat)
            if os.path.exists(p):
                os.remove(p)

        # 删除 bam（可选）
        if not keep_bam:
            bam_dir = os.path.join(args.align_dir, "hisat2file", srr)
            if os.path.exists(bam_dir):
                shutil.rmtree(bam_dir, ignore_errors=True)
                ts_log(f"[ChunkClean] {srr}: 已删除 bam 目录")


def _run_chunk(chunk_srrs, args):
    """
    对一个 chunk 的样本列表执行并行 process_sample，返回失败列表。
    """
    failed = []
    with ThreadPoolExecutor(max_workers=args.maxparallel) as executor:
        futures = {executor.submit(process_sample, srr, args): srr
                   for srr in chunk_srrs}
        for future, srr in futures.items():
            try:
                ok = future.result()
                if not ok:
                    failed.append(srr)
            except Exception as e:
                ts_log(f"[Error] {srr} 异常: {e}")
                failed.append(srr)
    return failed


def main():
    global HISAT2_INDEX, ANNOTATIONS

    args = parse_args()

    # 应用命令行覆盖：hisat2_index / anno_base
    if args.anno_base:
        ANNOTATIONS = _build_annotations(args.anno_base)
        ts_log(f"[Config] anno_base 已覆盖为: {args.anno_base}")
    if args.hisat2_index:
        HISAT2_INDEX[args.species] = args.hisat2_index
        ts_log(f"[Config] hisat2_index[{args.species}] 已覆盖为: {args.hisat2_index}")

    os.makedirs(args.fastq_dir, exist_ok=True)
    os.makedirs(args.qc_dir, exist_ok=True)
    os.makedirs(args.align_dir, exist_ok=True)

    # 获取该数据集所有 SRR
    sra_files = sorted(glob.glob(os.path.join(args.sra_dir, "*.sra")))
    if not sra_files:
        ts_log(f"[Warning] 在 {args.sra_dir} 中未找到任何 .sra 文件（该 GSE 可能无 RNA-seq 数据或下载失败）")
        # 写出 skipped marker 让 Snakemake 认为此任务已完成（跳过，不是真正失败）
        os.makedirs(os.path.dirname(args.output_marker), exist_ok=True)
        with open(args.output_marker, "w") as f:
            f.write(f"skipped: no .sra files found in {args.sra_dir}\n")
        ts_log(f"[Info] 已写出跳过标志文件: {args.output_marker}")
        sys.exit(0)

    srr_list = [os.path.splitext(os.path.basename(f))[0] for f in sra_files]
    chunk_size = args.sample_chunk_size

    ts_log(f"[Info] 数据集 {args.dataset_id} 共 {len(srr_list)} 个样本")
    ts_log(f"[Info] 并行数: {args.maxparallel}，每样本线程: {args.threads}")

    # ── 读取 cleanup 策略（从环境变量或默认值）──
    # Snakemake shell 中通过 export 传入（与原逻辑一致）
    keep_sra = os.environ.get("KEEP_SRA", "false").lower() == "true"
    keep_bam = os.environ.get("KEEP_BAM", "false").lower() == "true"

    all_failed = []

    if chunk_size <= 0:
        # ── 原始模式：一次性处理全部样本 ──
        ts_log(f"[Info] 不分块，一次性处理所有 {len(srr_list)} 个样本")
        ts_log(f"[Info] 样本列表: {srr_list}")
        all_failed = _run_chunk(srr_list, args)
    else:
        # ── 滚动窗口模式：每次处理 chunk_size 个样本 ──
        chunks = [srr_list[i:i + chunk_size]
                  for i in range(0, len(srr_list), chunk_size)]
        ts_log(f"[Info] 滚动窗口模式: chunk_size={chunk_size}，"
               f"共 {len(chunks)} 批，每批最多 {chunk_size} 个样本")

        for idx, chunk in enumerate(chunks, 1):
            ts_log(f"[Chunk {idx}/{len(chunks)}] 开始处理: {chunk}")
            failed = _run_chunk(chunk, args)

            # 只清理成功的样本，失败的保留 .sra 以供 --rerun-incomplete 重试
            success_srrs = [s for s in chunk if s not in failed]

            if failed:
                ts_log(f"[Chunk {idx}/{len(chunks)}] 以下样本失败（.sra 已保留供重试）: {failed}")
                all_failed.extend(failed)
            else:
                ts_log(f"[Chunk {idx}/{len(chunks)}] 全部成功")

            # 立即清理本批次中成功的样本
            if success_srrs:
                ts_log(f"[Chunk {idx}/{len(chunks)}] 清理成功样本磁盘占用: {success_srrs}")
                cleanup_chunk(success_srrs, args, keep_bam=keep_bam, keep_sra=keep_sra)
            ts_log(f"[Chunk {idx}/{len(chunks)}] 完成，继续下一批")

    # ── 最终结果 ──
    marker_dir = os.path.dirname(args.output_marker) or "."
    os.makedirs(marker_dir, exist_ok=True)
    failed_report = os.path.join(marker_dir, "failed_samples.txt")

    if all_failed:
        ts_log(f"[Error] 数据集 {args.dataset_id} 以下样本处理失败: {all_failed}")
        ts_log(f"[Error] 流程终止，请检查日志后重新运行（--rerun-incomplete）。")
        # 写出失败样本列表，方便 UI 展示和用户排查
        with open(failed_report, "w") as f:
            f.write(f"# {args.dataset_id} failed samples ({len(all_failed)}/{len(srr_list)})\n")
            f.write(f"# Generated at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            for srr in all_failed:
                f.write(f"{srr}\n")
        ts_log(f"[Info] 失败样本列表已写出: {failed_report}")
        # 不写 marker，以非零退出码通知 Snakemake 此任务失败
        sys.exit(1)
    else:
        ts_log(f"[Done] 数据集 {args.dataset_id} 所有样本处理完成")
        with open(args.output_marker, "w") as f:
            f.write("done\n")
        ts_log(f"[Info] 标志文件已写出: {args.output_marker}")
        # 清理旧的失败记录（如果存在）
        if os.path.exists(failed_report):
            os.remove(failed_report)
            ts_log(f"[Info] 已清除旧的 failed_samples.txt")


if __name__ == "__main__":
    main()
