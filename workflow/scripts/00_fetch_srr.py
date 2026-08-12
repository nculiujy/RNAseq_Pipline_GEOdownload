#!/usr/bin/env python3
"""
00_fetch_srr.py — 从 SRR_table.txt 自动爬取 SRR 信息

输入:
  workflow/resources/{species}/SRR_table.txt  — 每行一个 GSE 号（支持 # 注释）

输出（每个 GSE 独立子目录）:
  {outdir}/{GSE}/SraRunInfo.csv               — Run 详情
  {outdir}/{GSE}/SRR_Acc_List_rnaseq.txt      — RNA-seq SRR 列表
  {outdir}/{GSE}/SRR_Acc_List_all.txt         — 全部 SRR 列表
  {outdir}/GSE_SRR_summary.csv                — 汇总表（顶层）
  {outdir}/ALL_rnaseq_SRR.txt                 — 合并所有 GSE 的 RNA-seq SRR（顶层）

特性:
  - 幂等：已有 {GSE}/SraRunInfo.csv 的 GSE 自动跳过（--force 可强制重拉）
  - 复用 00_gse_to_srr.py 的核心爬取逻辑（use_subdir=True 模式）
  - 每个 GSE 独立目录，文件脉络清晰，无 GSE 前缀冗余

用法:
  # 从 SRR_table.txt 批量处理
  python workflow/scripts/00_fetch_srr.py \\
      --table  workflow/resources/homo/SRR_table.txt \\
      --outdir workflow/resources/homo/ \\
      [--force]

  # 单个 GSE（调试）
  python workflow/scripts/00_fetch_srr.py --gse GSE242225 \\
      --outdir workflow/resources/homo/
"""

import sys
import os
import re
import argparse


# ── 动态导入 00_gse_to_srr.py 的核心函数 ─────────────────────
def _import_gse_module():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    gse_script = os.path.join(script_dir, "00_gse_to_srr.py")
    import importlib.util
    spec = importlib.util.spec_from_file_location("gse_to_srr", gse_script)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_srr_table(table_path):
    """读取 SRR_table.txt，返回 GSE 号列表（去重保序，忽略注释行和空行）"""
    if not os.path.exists(table_path):
        print(f"[ERROR] SRR_table.txt 不存在: {table_path}")
        return []
    gses = []
    with open(table_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = re.search(r"(GSE\d+)", line, re.I)
            if m:
                gse = m.group(1).upper()
                if gse not in gses:
                    gses.append(gse)
    return gses


def parse_args():
    p = argparse.ArgumentParser(description="从 SRR_table.txt 自动爬取 SRR 信息（按 GSE 子目录）")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--table", help="SRR_table.txt 路径（每行一个 GSE 号）")
    group.add_argument("--gse",   help="单个 GSE 号（调试模式）")
    p.add_argument("--outdir", required=True,
                   help="输出根目录（每个 GSE 在此下创建子目录）")
    p.add_argument("--force",  action="store_true",
                   help="强制重新拉取（跳过缓存检查）")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 导入核心模块
    try:
        gse_mod = _import_gse_module()
    except Exception as e:
        print(f"[ERROR] 无法加载 00_gse_to_srr.py: {e}")
        sys.exit(1)

    # 确定 GSE 列表
    if args.gse:
        gses = [args.gse.strip().upper()]
        print(f"单 GSE 模式: {gses[0]}")
    else:
        gses = parse_srr_table(args.table)
        if not gses:
            print(f"[ERROR] SRR_table.txt 为空或无有效 GSE 号: {args.table}")
            sys.exit(1)
        print(f"从 {args.table} 读取 {len(gses)} 个 GSE")

    ok_count, skip_count, fail_list = 0, 0, []

    for gse in gses:
        gse_dir   = os.path.join(args.outdir, gse)
        info_path = os.path.join(gse_dir, "SraRunInfo.csv")

        # 幂等跳过
        if not args.force and os.path.exists(info_path) and os.path.getsize(info_path) > 0:
            skip_count += 1
            print(f"[Skip] {gse} → {gse_dir}/ （已有结果，--force 可重拉）")
            continue

        print(f"\n{'='*40}")
        print(f"[Fetch] {gse} → {gse_dir}/")
        print(f"{'='*40}")
        try:
            # use_subdir=True：写入 {outdir}/{GSE}/ 子目录，文件名不含 GSE 前缀
            ok, n_rna, msg = gse_mod.process_gse(gse, args.outdir, use_subdir=True)
            if ok:
                ok_count += 1
                print(f"[OK] {gse}: RNA-seq {n_rna} 个 SRR")
            else:
                fail_list.append((gse, msg))
                print(f"[FAIL] {gse}: {msg}")
        except Exception as e:
            fail_list.append((gse, str(e)))
            print(f"[ERROR] {gse}: {e}")

    print(f"\n{'='*40}")
    print(f"[汇总] 新拉取: {ok_count}  跳过: {skip_count}  失败: {len(fail_list)}")
    if fail_list:
        print("[失败列表]:")
        for g, m in fail_list:
            print(f"  失败: {g}: {m}")

    print(f"[Done] 输出目录: {args.outdir}")
    print(f"[Info] 每个 GSE 的数据已存放于 {args.outdir}/<GSE>/ 子目录下")


if __name__ == "__main__":
    main()
