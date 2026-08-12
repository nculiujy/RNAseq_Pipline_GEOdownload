#!/usr/bin/env python3
"""
05_merge_all_batches.py — 跨批次矩阵整合

扫描 result/{project}/*/03_Align_Filter/{species}/ 下所有 run_id 的 gene_abund.tab，
合并为全量表达矩阵。

输出:
  result/{project}/00_final_matrices/Matrices_*/*_matrix.csv
  result/{project}/00_final_matrices/sample_manifest.csv

用法:
  python 05_merge_all_batches.py \
      --result_base result/GBM_homo \
      --species homo \
      --output_dir result/GBM_homo/00_final_matrices
"""

import os
import sys
import argparse
import glob
import csv
from datetime import datetime
from collections import defaultdict


def ts_log(msg: str):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def parse_args():
    parser = argparse.ArgumentParser(description="Cross-batch matrix merge")
    parser.add_argument("--result_base", required=True,
                        help="result/{project} base directory")
    parser.add_argument("--species", required=True,
                        help="Species: homo or mouse")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for merged matrices")
    parser.add_argument("--cutoff", type=float, default=70.0,
                        help="Alignment rate cutoff for filtering (default: 70.0)")
    return parser.parse_args()


def scan_gene_abund_files(result_base, species):
    """
    递归扫描 result_base/*/03_Align_Filter/{species}/GSE*/注释类型/SRR/gene_abund.tab
    返回 dict: { annotation_type: [ (srr, run_id, gse, filepath) ] }
    """
    # 8 类注释子目录模式
    anno_patterns = [
        "mRNA/genecode/stringtie",
        "eRNA/EnhancerAtlas/stringtie",
        "eRNA/Ensembl/stringtie",
        "eRNA/FANTOM5/stringtie",
        "lncRNA/GENCODE/stringtie",
        "lncRNA/NONCODE/stringtie",
        "miRNA/miRBase/stringtie",
        "miRNA/MirGeneDB/stringtie",
    ]

    results = defaultdict(list)

    # 扫描所有 run_id 目录（排除 00_ 开头的项目级目录）
    run_id_dirs = sorted([
        d for d in glob.glob(os.path.join(result_base, "*"))
        if os.path.isdir(d) and not os.path.basename(d).startswith("00_")
    ])

    for run_dir in run_id_dirs:
        run_id = os.path.basename(run_dir)
        align_base = os.path.join(run_dir, "03_Align_Filter", species)
        if not os.path.isdir(align_base):
            continue

        # 扫描 GSE 目录
        gse_dirs = [d for d in glob.glob(os.path.join(align_base, "GSE*"))
                    if os.path.isdir(d)]

        for gse_dir in gse_dirs:
            gse = os.path.basename(gse_dir)

            for anno_path in anno_patterns:
                anno_dir = os.path.join(gse_dir, anno_path)
                if not os.path.isdir(anno_dir):
                    continue

                # 扫描 SRR 目录
                for srr_dir in glob.glob(os.path.join(anno_dir, "SRR*")):
                    if not os.path.isdir(srr_dir):
                        continue
                    srr = os.path.basename(srr_dir)
                    abund_file = os.path.join(srr_dir, "gene_abund.tab")
                    if os.path.exists(abund_file):
                        results[anno_path].append((srr, run_id, gse, abund_file))

    return results


def load_alignment_quality(result_base):
    """
    加载所有批次的 alignment_quality.csv，合并为 {sample_id: {rate, passed}} dict。
    """
    qc_data = {}
    for csv_path in glob.glob(os.path.join(result_base, "*", "03_Align_Filter", "alignment_quality.csv")):
        run_id = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(csv_path))))
        try:
            with open(csv_path, newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    sid = row.get("Sample_ID", "")
                    if sid:
                        qc_data[sid] = {
                            "rate": float(row.get("Alignment_Rate", 0)),
                            "passed": row.get("Passed", ""),
                            "run_id": run_id,
                        }
        except Exception as e:
            ts_log(f"[Warning] 读取 QC CSV 失败: {csv_path}: {e}")
    return qc_data


def merge_gene_abund(file_list, output_path):
    """
    合并多个 gene_abund.tab 文件为矩阵（gene × sample）。
    gene_abund.tab 格式: tab-separated, 列包含 Gene ID, Gene Name, Coverage, FPKM, TPM
    """
    import pandas as pd

    all_data = {}
    for srr, run_id, gse, filepath in file_list:
        try:
            df = pd.read_csv(filepath, sep="\t")
            # StringTie gene_abund.tab 列: Gene ID, Gene Name, Reference, Strand, Start, End, Coverage, FPKM, TPM
            if "Gene ID" in df.columns and "TPM" in df.columns:
                series = df.set_index("Gene ID")["TPM"]
                # 去重：同 SRR 出现在多个批次时保留最新（file_list 已按 run_id 排序）
                if srr not in all_data:
                    all_data[srr] = series
                else:
                    ts_log(f"[Warning] 样本 {srr} 在多个批次中出现，保留后者 (run_id={run_id})")
                    all_data[srr] = series
        except Exception as e:
            ts_log(f"[Warning] 读取失败: {filepath}: {e}")

    if not all_data:
        return None

    matrix = pd.DataFrame(all_data)
    matrix.index.name = "Gene_ID"
    matrix.fillna(0, inplace=True)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    matrix.to_csv(output_path)
    return matrix


def main():
    args = parse_args()

    ts_log(f"[Start] 跨批次合并 — result_base={args.result_base}, species={args.species}")

    # 1. 扫描所有 gene_abund.tab
    anno_files = scan_gene_abund_files(args.result_base, args.species)
    total_files = sum(len(v) for v in anno_files.values())
    ts_log(f"[Info] 扫描到 {total_files} 个 gene_abund.tab（{len(anno_files)} 类注释）")

    if total_files == 0:
        ts_log("[Warning] 未扫描到任何 gene_abund.tab，退出")
        sys.exit(0)

    # 2. 加载 QC 数据
    qc_data = load_alignment_quality(args.result_base)
    ts_log(f"[Info] 加载 {len(qc_data)} 条 QC 记录")

    # 3. 构建 sample_manifest
    manifest_rows = []
    all_srrs_seen = set()

    for anno_path, file_list in anno_files.items():
        for srr, run_id, gse, filepath in file_list:
            if srr in all_srrs_seen:
                continue
            all_srrs_seen.add(srr)

            qc_info = qc_data.get(srr, {})
            manifest_rows.append({
                "GSE": gse,
                "Run": srr,
                "Source_Batch": run_id,
                "Alignment_Rate": qc_info.get("rate", ""),
                "Passed": qc_info.get("passed", ""),
            })

    # 写 sample_manifest.csv
    os.makedirs(args.output_dir, exist_ok=True)
    manifest_path = os.path.join(args.output_dir, "sample_manifest.csv")
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["GSE", "Run", "Source_Batch", "Alignment_Rate", "Passed"])
        writer.writeheader()
        writer.writerows(manifest_rows)
    ts_log(f"[Info] sample_manifest.csv 已写出: {len(manifest_rows)} 条 → {manifest_path}")

    # 4. 按注释类型合并矩阵
    try:
        import pandas as pd
    except ImportError:
        ts_log("[Error] 需要 pandas: pip install pandas")
        sys.exit(1)

    for anno_path, file_list in anno_files.items():
        # 注释类型名（如 mRNA_genecode）
        anno_name = anno_path.replace("/stringtie", "").replace("/", "_")
        out_dir = os.path.join(args.output_dir, f"Matrices_all")
        out_file = os.path.join(out_dir, f"{anno_name}_matrix.csv")

        ts_log(f"[Merge] {anno_name}: {len(file_list)} 个样本")
        matrix = merge_gene_abund(file_list, out_file)
        if matrix is not None:
            ts_log(f"[Merge-Done] {anno_name}: {matrix.shape[0]} genes × {matrix.shape[1]} samples → {out_file}")
        else:
            ts_log(f"[Merge-Warning] {anno_name}: 无有效数据")

    ts_log("[Done] 跨批次合并完成")


if __name__ == "__main__":
    main()
