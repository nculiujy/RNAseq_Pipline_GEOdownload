"""
streamlit_app/core/state.py — pipeline_state.json 读写
"""
import os
import json
from datetime import datetime
import glob
import csv
import re


def load_pipeline_state(project, result_base="result"):
    """
    加载 result/{project}/pipeline_state.json，
    若不存在则扫描文件系统动态生成（慢，仅首次）
    """
    state_file = os.path.join(result_base, project, "pipeline_state.json")
    if os.path.exists(state_file):
        with open(state_file) as f:
            return json.load(f)
    return _build_state_from_fs(project, result_base)


def _build_state_from_fs(project, result_base):
    """扫描文件系统生成 pipeline_state（慢路径，首次或无缓存时使用）"""
    gses_info = []
    align_base = os.path.join(result_base, project, "03_Align_Filter")
    # 扫描所有 dataset_finished.txt
    for fpath in sorted(glob.glob(
            os.path.join(align_base, "**", "dataset_finished.txt"), recursive=True)):
        # 路径结构: 03_Align_Filter/{species}/{gse}/dataset_finished.txt
        rel = os.path.relpath(fpath, align_base)
        parts = rel.split(os.sep)
        if len(parts) >= 2:
            gse = parts[-2]  # 倒数第二级是 GSE
        else:
            continue
        with open(fpath) as f:
            content = f.read()
        content_lower = content.lower()
        if "failed" in content_lower:
            status = "failed"
            failed = re.findall(r"'(SRR\d+)'", content)
        elif content_lower.startswith("skipped"):
            status = "skipped"
            failed = []
        else:
            status = "done"
            failed = []
        gses_info.append({
            "gse": gse,
            "status": status,
            "failed_samples": failed,
            "marker": fpath
        })

    # 读取 alignment_quality.csv（若存在）
    qc_csv = os.path.join(align_base, "alignment_quality.csv")
    align_rates = {}
    if os.path.exists(qc_csv):
        with open(qc_csv, newline="") as f:
            for row in csv.DictReader(f):
                srr = row.get("Sample_ID", "")
                rate = float(row.get("Alignment_Rate", 0) or 0)
                align_rates[srr] = rate

    # 补充 align_rate_mean 到每个 GSE
    for g in gses_info:
        gse = g["gse"]
        gse_rates = [v for k, v in align_rates.items() if k.startswith("SRR")]
        g["align_rate_mean"] = round(sum(gse_rates) / len(gse_rates), 1) if gse_rates else None

    # 磁盘信息
    import shutil
    disk = shutil.disk_usage(result_base)
    disk_info = {
        "free_gb":  round(disk.free / 1e9, 1),
        "total_gb": round(disk.total / 1e9, 1),
        "used_gb":  round(disk.used / 1e9, 1),
    }

    state = {
        "generated_at": datetime.now().isoformat(),
        "projects": [{
            "name":  project,
            "gses":  gses_info,
            "disk":  disk_info,
        }]
    }
    return state


def save_pipeline_state(state, project, result_base="result"):
    """将 state 写回 pipeline_state.json"""
    os.makedirs(os.path.join(result_base, project), exist_ok=True)
    state_file = os.path.join(result_base, project, "pipeline_state.json")
    with open(state_file, "w") as f:
        json.dump(state, f, ensure_ascii=False, indent=2)


def get_project_state(state, project):
    """从完整 state 中提取指定 project 的数据"""
    for proj in state.get("projects", []):
        if proj.get("name") == project:
            return proj
    return {}


def list_projects(result_base="result"):
    """扫描 result/ 下所有项目目录"""
    if not os.path.isdir(result_base):
        return []
    return [d for d in sorted(os.listdir(result_base))
            if os.path.isdir(os.path.join(result_base, d))
            and not d.startswith(".")]
