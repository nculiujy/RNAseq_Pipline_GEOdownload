"""
RNAseq_GEO — Snakemake 主入口
==============================
流程概述（01-04，GSE 列表通过 UI 输入并用 00_fetch_srr.py 爬取）:
  01. download_gse     per-GSE prefetch 下载 .sra（resources: gse_slots=1 串行）
  02. dataset_pipeline 每个 GSE 数据集: fasterq-dump → fastp → HISAT2 → StringTie
  03. filter_align     扫描比对率日志 → alignment_quality.csv
  04. merge_matrices   合并所有 GSE 表达量矩阵 → 最终表达矩阵

用法:
  # 全量跑（读 SRR_table.txt 全部 32 个 GSE）
  bash run_all.sh

  # 批次跑（用指定 txt 文件，结果存到带时间戳的子目录）
  snakemake -j 32 --resources gse_slots=1 \\
      --config batch_file=workflow/resources/homo/batch01_5gse.txt \\
               run_id=20260808_batch01

  # dry-run
  snakemake -n
"""

import os
import re
from datetime import datetime

configfile: "config/config.yaml"

# ──────────────────────────────────────────────
# run_id：结果子目录名 = {pname}/{run_id}/
# 优先级: --config run_id=xxx > 自动生成（时间戳+批次文件名）
# ──────────────────────────────────────────────
_batch_file = config.get("batch_file", None)   # 指定批次文件
_run_id     = config.get("run_id",     None)   # 指定结果子目录名

if not _run_id:
    if _batch_file:
        _batch_name = os.path.splitext(os.path.basename(_batch_file))[0]
        _run_id = f"{datetime.now().strftime('%Y%m%d')}_{_batch_name}"
    else:
        _run_id = datetime.now().strftime('%Y%m%d_all')

print(f"[Snakefile] run_id = {_run_id}")

# ──────────────────────────────────────────────
# 读取 GSE 列表
# 优先级: --config batch_file=xxx > batch_input.txt > SRR_table.txt
# ──────────────────────────────────────────────

def _load_gse_list(rawdata_dir, species):
    """读取 GSE 列表（返回 GSE 号列表）"""
    # 优先级 1: --config batch_file 指定文件
    if _batch_file and os.path.exists(_batch_file):
        input_path = _batch_file
        src = f"--config batch_file"
    # 优先级 2: batch_input.txt（UI 旧版兼容）
    elif os.path.exists(os.path.join(rawdata_dir, "batch_input.txt")):
        input_path = os.path.join(rawdata_dir, "batch_input.txt")
        src = "batch_input.txt"
    # 优先级 3: SRR_table.txt（全量）
    else:
        input_path = os.path.join(rawdata_dir, "SRR_table.txt")
        src = "SRR_table.txt"

    gse_list = []
    if not os.path.exists(input_path):
        print(f"[Snakefile] WARNING: 未找到 GSE 列表文件: {input_path}")
        return gse_list

    with open(input_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = re.search(r"(GSE\d+)", line, re.I)
            if m:
                gse = m.group(1).upper()
                if gse not in gse_list:
                    gse_list.append(gse)

    print(f"[Snakefile] GSE 列表来源: {input_path} ({src}, {len(gse_list)} 个 GSE)")
    return gse_list


# 将 run_id 注入 config（让 smk rules 可以通过 config["run_id"] 读取）
config["run_id"] = _run_id

# ──────────────────────────────────────────────
# 包含各模块规则文件
# ──────────────────────────────────────────────
include: "workflow/rules/01_download_sra.smk"
include: "workflow/rules/02_dataset_pipeline.smk"
include: "workflow/rules/03_filter_alignment.smk"
include: "workflow/rules/04_merge_matrices.smk"
include: "workflow/rules/05_merge_all_batches.smk"

# ──────────────────────────────────────────────
# 磁盘告警钩子
# ──────────────────────────────────────────────
import shutil as _shutil

onstart:
    _min_free = config.get("batch", {}).get("min_free_gb", 300)
    _free_gb  = _shutil.disk_usage(".").free / 1e9
    if _free_gb < _min_free:
        print(f"\n{'!'*50}")
        print(f"⚠️  磁盘告警: 剩余 {_free_gb:.1f} GB < 阈值 {_min_free} GB")
        print(f"{'!'*50}\n")
        # 严格模式：磁盘不足时中止（防止写满磁盘导致数据损坏）
        _strict_disk = config.get("batch", {}).get("strict_disk_check", False)
        if _strict_disk:
            raise RuntimeError(
                f"磁盘空间不足: {_free_gb:.1f} GB < {_min_free} GB。"
                f"请释放磁盘空间或在 config.yaml batch.min_free_gb 中降低阈值。"
                f"若要强制启动，设置 batch.strict_disk_check: false"
            )
    # 并发预检
    _pt = config.get("pipeline_threads", 8)
    _pp = config.get("pipeline_parallel", 4)
    _gs = config.get("batch", {}).get("gse_slots", 1)
    _cpu = os.cpu_count() or 1
    _total_needed = _pt * _pp * _gs
    if _total_needed > _cpu:
        print(f"[onstart] ⚠️  并发偏高: {_gs}×{_pp}×{_pt}={_total_needed} > CPU核心数 {_cpu}")
    print(f"[onstart] run_id = {_run_id}, disk_free = {_free_gb:.1f} GB")

# ──────────────────────────────────────────────
# 构建最终目标文件列表
# 结果目录：result/{pname}/{run_id}/
# ──────────────────────────────────────────────
TARGET_FILES = []

PROJECTS = config.get("projects", [])

for proj in PROJECTS:
    pname      = proj["project_name"]
    species    = proj["species"]
    rawdata    = proj.get("rawdata_dir", f"workflow/resources/{species}")
    modules    = proj.get("modules", {})

    # ★ 将 run_id 合并入 project wildcard（无需修改任何 smk 规则文件）
    # 结果路径：result/GBM_homo/20260808_batch01/03_Align_Filter/...
    # {project} wildcard 在规则中匹配 "GBM_homo/20260808_batch01"
    effective_project = f"{pname}/{_run_id}"

    gse_list = _load_gse_list(rawdata, species)

    # 01+02 per-GSE（DAG 自动推导 01→02）
    if modules.get("02_dataset_pipeline", False):
        for gse in gse_list:
            TARGET_FILES.append(
                os.path.join("result", effective_project, "03_Align_Filter", species, gse, "dataset_finished.txt")
            )

    # 03 比对质量过滤
    if modules.get("03_filter_alignment", False):
        TARGET_FILES.append(
            os.path.join("result", effective_project, "03_Align_Filter", "alignment_quality.csv")
        )

    # 04 合并矩阵
    if modules.get("04_merge_matrices", False):
        TARGET_FILES.append(
            os.path.join("result", effective_project, "04_merge_matrices", "Merge_finished.txt")
        )

rule all:
    input:
        TARGET_FILES
