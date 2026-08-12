"""
streamlit_app/core/geo.py — GEO/SRA 信息读取
"""
import os
import csv
import glob
import json


def load_all_sra_info(sra_info_dir):
    """
    扫描 sra_info_dir 下所有 SRR 信息文件，合并为含 GSE 列的大表。
    支持两种目录结构：
      - 新结构（按 GSE 子目录）: {sra_info_dir}/{GSE}/SraRunInfo.csv
      - 旧结构（平铺）: {sra_info_dir}/SraRunInfo_{GSE}.csv
    """
    rows = []
    # 新结构：{GSE}/SraRunInfo.csv
    for gse_dir in sorted(glob.glob(os.path.join(sra_info_dir, "GSE*"))):
        gse = os.path.basename(gse_dir)
        fpath = os.path.join(gse_dir, "SraRunInfo.csv")
        if not os.path.exists(fpath):
            continue
        try:
            with open(fpath, newline="") as f:
                for row in csv.DictReader(f):
                    row["GSE"] = gse
                    rows.append(row)
        except Exception as e:
            print(f"[WARN] 读取 {fpath} 失败: {e}")
    # 旧结构（向下兼容）: SraRunInfo_{GSE}.csv
    if not rows:
        for fpath in sorted(glob.glob(os.path.join(sra_info_dir, "SraRunInfo_*.csv"))):
            gse = os.path.basename(fpath).replace("SraRunInfo_", "").replace(".csv", "")
            try:
                with open(fpath, newline="") as f:
                    for row in csv.DictReader(f):
                        row["GSE"] = gse
                        rows.append(row)
            except Exception as e:
                print(f"[WARN] 读取 {fpath} 失败: {e}")
    return rows


def load_gse_summary(sra_info_dir):
    """读取 GSE_SRR_summary.csv，返回 list of dict"""
    summary_csv = os.path.join(sra_info_dir, "GSE_SRR_summary.csv")
    if not os.path.exists(summary_csv):
        return []
    with open(summary_csv, newline="") as f:
        return list(csv.DictReader(f))


def load_llm_card(project, gse, result_base="result"):
    """加载 LLM 解读卡 JSON，返回 dict 或 None"""
    path = os.path.join(result_base, project, "00_data_intel", f"{gse}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def load_llm_card_md(project, gse, result_base="result"):
    """加载 LLM 解读卡 Markdown，返回字符串或 None"""
    path = os.path.join(result_base, project, "00_data_intel", f"{gse}.md")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return f.read()
