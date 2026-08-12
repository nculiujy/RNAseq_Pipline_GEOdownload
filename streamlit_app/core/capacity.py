"""
streamlit_app/core/capacity.py — 容量规划数据读取
"""
import os
import json


def load_plan(project, result_base="result"):
    """读取 plan.json，返回 summary dict 或 None"""
    plan_file = os.path.join(result_base, project, "00_planning", "plan.json")
    if not os.path.exists(plan_file):
        return None
    with open(plan_file) as f:
        return json.load(f)


def load_plan_csv(project, result_base="result"):
    """读取 plan.csv，返回 pandas DataFrame 或 None"""
    try:
        import pandas as pd
        csv_file = os.path.join(result_base, project, "00_planning", "plan.csv")
        if not os.path.exists(csv_file):
            return None
        return pd.read_csv(csv_file)
    except ImportError:
        return None


def get_gantt_html_path(project, result_base="result"):
    """返回甘特图 HTML 路径（若存在）"""
    path = os.path.join(result_base, project, "00_planning", "gantt.html")
    return path if os.path.exists(path) else None
