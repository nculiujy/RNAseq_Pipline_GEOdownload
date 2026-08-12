#!/usr/bin/env python3
"""
10_capacity_planner.py — 容量规划与批次预估

功能:
  1. 读取 SraRunInfo_*.csv，汇总每个 GSE 的样本数、总 spots/bases
  2. 估算每 GSE 的磁盘需求（.sra / 中间 fastq / 最终产物）
  3. 按磁盘安全阈值确定并发/批次方案
  4. 预估总耗时（L0 先验 / L1 基准校准 / L2 滚动更新）
  5. 输出 plan.json + plan.csv + 可选甘特图 HTML

用法:
  # 基础规划（L0 先验）
  python workflow/scripts/10_capacity_planner.py \\
      --sra_info   workflow/resources/homo/ \\
      --config     config/config.yaml \\
      --output_dir result/GBM_homo/00_planning

  # 带 L1 基准校准（提供已跑完样本的日志目录）
  python workflow/scripts/10_capacity_planner.py \\
      --sra_info   workflow/resources/homo/ \\
      --config     config/config.yaml \\
      --bench_log  result/GBM_homo/03_Align_Filter \\
      --output_dir result/GBM_homo/00_planning
"""

import os
import sys
import csv
import json
import math
import glob
import argparse
import datetime
from collections import defaultdict

# ─────────────────────────────────────────────
# L0 先验系数（来自文档实测系数）
# ─────────────────────────────────────────────
L0_SRA_PER_BASE     = 0.87   # sra 文件 / bases（比例）
L0_FASTQ_PER_BASE   = 2.0    # 原始 fastq / bases
L0_CLEAN_PER_BASE   = 1.9    # clean fastq / bases
L0_BAM_PER_BASE     = 0.375  # dedup bam / bases（1.5G / 4G 约）
L0_GTF_PER_SAMPLE   = 50e6   # StringTie gtf 约 50 MB/样本（8 注释合计）

# 时间模型：T_样本(s) ≈ alpha × spots_M + beta
L0_ALPHA_S_PER_M   = 0.002 * 60  # 0.002 min/spots_M → 秒
L0_BETA_S          = 8.0 * 60    # 8 min 固定开销 → 秒
IO_OVERHEAD        = 0.15        # IO/网络损耗系数


def parse_args():
    p = argparse.ArgumentParser(description="容量规划与批次预估")
    p.add_argument("--sra_info",   required=True,
                   help="SraRunInfo_*.csv 所在目录（或含 GSE_SRR_summary.csv 的目录）")
    p.add_argument("--config",     default="config/config.yaml",
                   help="config.yaml 路径（读取 pipeline_parallel/threads/planning 字段）")
    p.add_argument("--bench_log",  default=None,
                   help="已完成的 03_Align_Filter 目录，用于 L1/L2 时间校准")
    p.add_argument("--output_dir", default="result/planning",
                   help="输出目录（plan.json / plan.csv / gantt.html）")
    p.add_argument("--project",    default="",
                   help="项目名（写入 plan.json，用于 UI 区分多项目）")
    return p.parse_args()


# ──────────────────────────────────────────────
# 工具函数
# ──────────────────────────────────────────────

def load_config(config_path):
    """简单读取 YAML（避免依赖 pyyaml，支持基础 key:value 格式）"""
    cfg = {}
    try:
        with open(config_path) as f:
            import yaml
            cfg = yaml.safe_load(f) or {}
    except ImportError:
        # fallback: 逐行读基础 key: value
        with open(config_path) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and ":" in line:
                    k, _, v = line.partition(":")
                    cfg[k.strip()] = v.strip()
    return cfg


def load_sra_info_files(sra_info_dir):
    """
    扫描 sra_info_dir 下所有 GSE 子目录中的 SraRunInfo.csv，
    返回 {gse: [{"Run": ..., "spots": ..., "bases": ..., ...}, ...]}

    目录结构: {sra_info_dir}/{GSE}/SraRunInfo.csv
    """
    gse_data = {}
    # 新结构：{GSE}/SraRunInfo.csv
    for gse_dir in sorted(glob.glob(os.path.join(sra_info_dir, "GSE*"))):
        if not os.path.isdir(gse_dir):
            continue
        gse = os.path.basename(gse_dir)
        fpath = os.path.join(gse_dir, "SraRunInfo.csv")
        if not os.path.exists(fpath):
            continue
        rows = []
        try:
            with open(fpath, newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    rows.append(row)
        except Exception as e:
            print(f"[WARN] 读取 {fpath} 失败: {e}")
        if rows:
            gse_data[gse] = rows
    return gse_data


def load_bench_times(bench_dir):
    """
    扫描 bench_dir 下所有 QC_results.log，提取文件 mtime 作为粗略耗时参考。
    未来：解析 [2026-08-08 10:00:01] 时间戳日志 → 精确各步耗时。
    返回 {srr: {"total_s": float}} 或空字典
    """
    times = {}
    for fpath in glob.glob(os.path.join(bench_dir, "**", "QC_results.log"), recursive=True):
        srr = os.path.basename(os.path.dirname(fpath))
        bam = fpath.replace("QC_results.log", f"{srr}.dedup.bam")
        if os.path.exists(bam):
            # 用 dedup bam 和 QC log 的 mtime 差估算
            t_start = os.path.getmtime(fpath)
            t_end   = os.path.getmtime(bam)
            if t_end > t_start:
                times[srr] = {"total_s": t_end - t_start}
    return times


def calibrate_alpha_beta(bench_times, gse_data_flat):
    """
    L1 校准：用实测数据拟合 T = alpha * spots_M + beta
    返回 (alpha, beta) 或 L0 默认值
    """
    X, Y = [], []
    srr_spots = {r["Run"]: int(r.get("spots") or 0)
                 for rows in gse_data_flat.values()
                 for r in rows if r.get("Run")}
    for srr, t in bench_times.items():
        spots = srr_spots.get(srr, 0)
        if spots > 0 and t["total_s"] > 0:
            X.append(spots / 1e6)
            Y.append(t["total_s"])
    if len(X) < 3:
        return L0_ALPHA_S_PER_M, L0_BETA_S  # 样本不足，用先验
    # 简单最小二乘
    n   = len(X)
    sx  = sum(X)
    sy  = sum(Y)
    sxx = sum(xi ** 2 for xi in X)
    sxy = sum(X[i] * Y[i] for i in range(n))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-9:
        return L0_ALPHA_S_PER_M, L0_BETA_S
    alpha = (n * sxy - sx * sy) / denom
    beta  = (sy - alpha * sx) / n
    return max(alpha, 0), max(beta, 0)


def estimate_gse(gse, rows, alpha, beta, parallel, disk_free_bytes, safety_ratio):
    """
    计算单个 GSE 的估算指标
    返回 dict
    """
    n_samples  = len(rows)
    total_spots = sum(int(r.get("spots") or 0) for r in rows)
    total_bases = sum(int(r.get("bases") or 0) for r in rows)
    avg_spots   = total_spots / n_samples if n_samples else 0
    avg_bases   = total_bases / n_samples if n_samples else 0

    # 存储估算（bytes）
    sra_size     = total_bases * L0_SRA_PER_BASE
    final_size   = sra_size + n_samples * L0_GTF_PER_SAMPLE  # 保留产物
    peak_1sample = avg_bases * (L0_SRA_PER_BASE + L0_FASTQ_PER_BASE
                                + L0_CLEAN_PER_BASE + L0_BAM_PER_BASE)
    peak_size    = min(parallel, n_samples) * peak_1sample

    # 时间估算
    t_sample_s  = alpha * (avg_spots / 1e6) + beta
    n_rounds    = math.ceil(n_samples / parallel) if parallel > 0 else n_samples
    t_gse_s     = n_rounds * t_sample_s * (1 + IO_OVERHEAD)

    # 磁盘是否足够（峰值 < 安全阈值）
    disk_ok = peak_size < disk_free_bytes * safety_ratio

    return {
        "gse":         gse,
        "n_samples":   n_samples,
        "total_spots": total_spots,
        "total_bases": total_bases,
        "sra_gb":      round(sra_size / 1e9, 2),
        "final_gb":    round(final_size / 1e9, 2),
        "peak_gb":     round(peak_size / 1e9, 2),
        "t_gse_h":     round(t_gse_s / 3600, 2),
        "disk_ok":     disk_ok,
    }


def get_disk_free(path="."):
    """返回 path 所在分区的可用字节数"""
    import shutil
    usage = shutil.disk_usage(path)
    return usage.free


def load_pipeline_state(project, result_base="result"):
    """读取 pipeline_state.json（若存在），返回每个 GSE 的运行状态"""
    state_file = os.path.join(result_base, project, "pipeline_state.json")
    if not os.path.exists(state_file):
        return {}
    with open(state_file) as f:
        state = json.load(f)
    gse_status = {}
    for proj in state.get("projects", []):
        if proj.get("name") == project:
            for g in proj.get("gses", []):
                gse_status[g["gse"]] = g.get("status", "unknown")
    return gse_status


def build_gantt_html(plan_rows, output_path):
    """
    生成简单的 HTML 甘特图（使用内联 plotly CDN，不需要本地 plotly 安装）
    """
    now = datetime.datetime.now()
    tasks, starts, durations, statuses = [], [], [], []
    cursor = now
    for row in plan_rows:
        tasks.append(row["gse"])
        starts.append(cursor.isoformat())
        durations.append(row["t_gse_h"])
        statuses.append(row.get("status", "pending"))
        cursor += datetime.timedelta(hours=row["t_gse_h"])

    # 颜色映射
    color_map = {"done": "#4CAF50", "running": "#2196F3",
                 "failed": "#f44336", "pending": "#9E9E9E"}
    colors = [color_map.get(s, "#9E9E9E") for s in statuses]

    # 生成简单 HTML 表格式甘特图（不依赖 plotly）
    rows_html = ""
    for i, row in enumerate(plan_rows):
        color = colors[i]
        status_emoji = {"done": "✅", "running": "⏳", "failed": "❌", "pending": "⬜"}.get(
            statuses[i], "⬜")
        rows_html += f"""
        <tr>
          <td>{row['gse']}</td>
          <td>{row['n_samples']}</td>
          <td>{row['sra_gb']:.1f} GB</td>
          <td>{row['final_gb']:.1f} GB</td>
          <td>{row['peak_gb']:.1f} GB</td>
          <td>{row['t_gse_h']:.1f} h</td>
          <td><span style="color:{color}">{status_emoji} {statuses[i]}</span></td>
          <td>{'✅' if row['disk_ok'] else '⚠️'}</td>
        </tr>"""

    total_h = sum(r["t_gse_h"] for r in plan_rows)
    total_samples = sum(r["n_samples"] for r in plan_rows)
    total_sra_gb = sum(r["sra_gb"] for r in plan_rows)
    total_final_gb = sum(r["final_gb"] for r in plan_rows)

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<title>RNAseq_GEO 批次规划</title>
<style>
  body {{ font-family: sans-serif; margin: 20px; }}
  h1 {{ color: #333; }}
  .kpi {{ display: flex; gap: 20px; margin-bottom: 20px; flex-wrap: wrap; }}
  .kpi-card {{ background: #f5f5f5; border-radius: 8px; padding: 16px; min-width: 150px; }}
  .kpi-card .val {{ font-size: 1.8em; font-weight: bold; color: #1976D2; }}
  .kpi-card .label {{ color: #666; font-size: 0.9em; }}
  table {{ border-collapse: collapse; width: 100%; }}
  th {{ background: #1976D2; color: white; padding: 8px 12px; text-align: left; }}
  td {{ padding: 6px 12px; border-bottom: 1px solid #eee; }}
  tr:hover td {{ background: #f9f9f9; }}
  .warn {{ color: #e65100; }}
  .generated {{ color: #999; font-size: 0.85em; margin-top: 20px; }}
</style>
</head>
<body>
<h1>📊 RNAseq_GEO 批次规划报告</h1>
<div class="kpi">
  <div class="kpi-card"><div class="val">{total_samples}</div><div class="label">总样本数</div></div>
  <div class="kpi-card"><div class="val">{total_sra_gb:.0f} GB</div><div class="label">总 SRA 大小（估）</div></div>
  <div class="kpi-card"><div class="val">{total_final_gb:.0f} GB</div><div class="label">最终产物（估）</div></div>
  <div class="kpi-card"><div class="val">{total_h:.1f} h</div><div class="label">预估总耗时</div></div>
  <div class="kpi-card"><div class="val">{len(plan_rows)}</div><div class="label">GSE 批次数</div></div>
</div>
<table>
  <thead>
    <tr>
      <th>GSE</th><th>样本数</th><th>SRA(GB)</th><th>产物(GB)</th>
      <th>峰值(GB)</th><th>预估耗时</th><th>状态</th><th>磁盘OK</th>
    </tr>
  </thead>
  <tbody>
    {rows_html}
  </tbody>
</table>
<p class="generated">生成时间: {now.strftime('%Y-%m-%d %H:%M:%S')} | 精度: L0 先验 ±50%（首批完成后自动校准）</p>
</body>
</html>"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[Plan] 甘特图已写出: {output_path}")


def main():
    args = parse_args()

    # ── 读取配置 ──
    cfg = load_config(args.config) if os.path.exists(args.config) else {}
    planning = cfg.get("planning", {}) if isinstance(cfg.get("planning"), dict) else {}
    parallel = int(cfg.get("pipeline_parallel", planning.get("parallel", 4)))
    safety_ratio = float(planning.get("safety_disk_ratio", 0.8))
    alpha_l0 = float(planning.get("alpha_per_m_spots", L0_ALPHA_S_PER_M))
    beta_l0  = float(planning.get("beta_fixed_min", L0_BETA_S / 60)) * 60  # 转秒

    # ── 读取 SraRunInfo ──
    gse_data = load_sra_info_files(args.sra_info)
    if not gse_data:
        print(f"[ERROR] 在 {args.sra_info} 中未找到 SraRunInfo_*.csv 文件")
        sys.exit(1)
    print(f"[Plan] 加载 {len(gse_data)} 个 GSE 的 SRA 信息")

    # ── L1/L2 校准（可选）──
    alpha, beta = alpha_l0, beta_l0
    if args.bench_log and os.path.isdir(args.bench_log):
        bench_times = load_bench_times(args.bench_log)
        if bench_times:
            alpha, beta = calibrate_alpha_beta(bench_times, gse_data)
            print(f"[Plan] L1 校准: α={alpha:.4f} s/spots_M, β={beta:.1f} s "
                  f"（基于 {len(bench_times)} 个已完成样本）")
        else:
            print("[Plan] bench_log 中未找到已完成样本，使用 L0 先验")

    # ── 磁盘剩余 ──
    disk_free = get_disk_free(args.output_dir if os.path.exists(args.output_dir) else ".")
    disk_free_gb = disk_free / 1e9

    # ── 读取已有 pipeline_state（更新运行状态）──
    gse_status = load_pipeline_state(args.project) if args.project else {}

    # ── 计算每 GSE 的估算 ──
    plan_rows = []
    for gse, rows in sorted(gse_data.items()):
        est = estimate_gse(gse, rows, alpha, beta, parallel, disk_free, safety_ratio)
        est["status"] = gse_status.get(gse, "pending")
        plan_rows.append(est)

    # 汇总
    total_h = sum(r["t_gse_h"] for r in plan_rows)
    total_samples = sum(r["n_samples"] for r in plan_rows)
    total_sra_gb  = sum(r["sra_gb"] for r in plan_rows)
    total_final_gb = sum(r["final_gb"] for r in plan_rows)

    # ── 输出 ──
    os.makedirs(args.output_dir, exist_ok=True)

    # plan.csv
    csv_path = os.path.join(args.output_dir, "plan.csv")
    fieldnames = ["gse", "n_samples", "total_spots", "total_bases",
                  "sra_gb", "final_gb", "peak_gb", "t_gse_h", "disk_ok", "status"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(plan_rows)
    print(f"[Plan] 计划表已写出: {csv_path}")

    # plan.json
    json_path = os.path.join(args.output_dir, "plan.json")
    summary = {
        "generated_at":   datetime.datetime.now().isoformat(),
        "project":        args.project,
        "n_gses":         len(plan_rows),
        "total_samples":  total_samples,
        "total_sra_gb":   round(total_sra_gb, 1),
        "total_final_gb": round(total_final_gb, 1),
        "disk_free_gb":   round(disk_free_gb, 1),
        "disk_usage_pct": round(total_final_gb / disk_free_gb * 100, 1) if disk_free_gb else 0,
        "total_hours":    round(total_h, 1),
        "alpha":          alpha,
        "beta":           beta,
        "calibration_level": "L1" if args.bench_log else "L0",
        "gses":           plan_rows
    }
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"[Plan] JSON 已写出: {json_path}")

    # gantt.html
    gantt_path = os.path.join(args.output_dir, "gantt.html")
    build_gantt_html(plan_rows, gantt_path)

    # 控制台摘要
    print(f"\n{'='*50}")
    print(f"[Plan] 摘要（项目: {args.project or 'N/A'}）")
    print(f"  GSE 批次数 : {len(plan_rows)}")
    print(f"  总样本数   : {total_samples}")
    print(f"  总 SRA 估算: {total_sra_gb:.1f} GB")
    print(f"  最终产物   : {total_final_gb:.1f} GB")
    print(f"  磁盘剩余   : {disk_free_gb:.1f} GB（占用 {total_final_gb/disk_free_gb*100:.0f}%）")
    print(f"  预估总耗时 : {total_h:.1f} h ≈ {total_h/24:.1f} 天（L0 先验，精度 ±50%）")
    print(f"{'='*50}")
    disk_warn = total_final_gb / disk_free_gb if disk_free_gb else 0
    if disk_warn > 0.5:
        print(f"  ⚠️  磁盘警告：最终产物占磁盘剩余的 {disk_warn*100:.0f}%，建议扩容或清理！")


if __name__ == "__main__":
    main()
