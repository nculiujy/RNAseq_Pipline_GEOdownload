#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gse_to_srr.py — Pipeline 第 0 步：GSE 列表 -> SRA Run 列表

链路: GSE --(GEO SOFT, FTP)--> BioProject --(E-utilities db=sra)--> SRR runs
输入: 任意含 GSE_ID 列的 CSV (如 GSE_list.csv / GSE_SRR_summary.csv)
输出: 每个 GSE 一份 SRR_Acc_List_<GSE>_rnaseq.txt + SraRunInfo_<GSE>.csv，
      以及汇总 GSE_SRR_summary.csv 和合并列表 ALL_rnaseq_SRR.txt
特性: 增量幂等 —— 已有结果的 GSE 自动跳过，重跑只需几秒；--force 强制重拉。
用法:
    python3 gse_to_srr.py --csv GSE_SRR_summary.csv [-o out/] [--force]
    python3 gse_to_srr.py GSE242225                    # 单个 GSE
依赖: Python3 标准库 (urllib/gzip/json)，无需第三方包
注意: 纯芯片数据(无 SRA)会跳过并记录；strategy 非 RNA 的 run 会被过滤到 _rnaseq 列表外。
"""

import sys, os, re, io, json, time, csv, gzip, glob
import urllib.request, urllib.parse
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP_SOFT = "https://ftp.ncbi.nlm.nih.gov/geo/series/{dir}/{gse}/soft/{gse}_family.soft.gz"
USER_AGENT = "gse_to_srr/1.0 (mailto:{})".format(
    os.environ.get("GEO_CONTACT_EMAIL", "your_email@example.com")
)
SLEEP = 0.4          # eutils 限速 ~3 req/s，请求间间隔(秒)
MAX_RETRY = 3        # 单请求最大重试次数

def fetch(url, retry=MAX_RETRY, post=None):
    for i in range(retry):
        try:
            req = urllib.request.Request(url, data=post, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception as e:
            if i == retry - 1:
                raise
            time.sleep(2 * (i + 1))
    raise RuntimeError("fetch failed: " + url)

def eutils(path, params, retmode="xml"):
    params = {k: v for k, v in params.items() if v is not None}
    qs = urllib.parse.urlencode(params)
    raw = fetch(f"{EUTILS}/{path}?{qs}&retmode={retmode}")
    time.sleep(SLEEP)
    return raw

def esearch(db, term, retmax=100000):
    xml = eutils("esearch.fcgi", {"db": db, "term": term, "retmax": retmax}).decode()
    ids = re.findall(r"<Id>(\d+)</Id>", xml)
    count = re.search(r"<Count>(\d+)</Count>", xml)
    err = re.search(r"<ERROR>(.*?)</ERROR>", xml)
    if err:
        raise RuntimeError(f"esearch error: {err.group(1)}")
    return ids, int(count.group(1)) if count else len(ids)

def soft_text(gse):
    """下载 GSE 的 family.soft.gz，返回解压文本。目录规则: GSE + (编号//1000) + 'nnn'"""
    m = re.fullmatch(r"GSE(\d+)", gse.upper())
    if not m:
        raise ValueError(f"不是合法的GSE号: {gse}")
    digits = m.group(1)
    # 正确规则: GSE46523 → GSE46nnn（整除1000）；旧版 digits[:3] 在 5 位 GSE 号上会算错
    url = FTP_SOFT.format(dir=f"GSE{int(digits) // 1000}nnn", gse=gse.upper())
    print(f"[1/4] 下载 SOFT: {url}")
    raw = fetch(url)
    if raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw).decode("utf-8", "replace")
    return raw.decode("utf-8", "replace")

def get_bioproject(text):
    """从 SOFT 提取 BioProject (Series_relation)。返回 PRJNAxxx 或 None"""
    for line in text.splitlines():
        if "BioProject:" in line:
            m = re.search(r"PRJN[A-Z]?\d+", line)
            if m:
                return m.group(0)
    # fallback: 任意含 PRJNA 的行
    m = re.search(r"PRJN[A-Z]?\d+", text)
    return m.group(0) if m else None

def parse_esummary(json_data):
    """解析 db=sra esummary JSON -> [{run, srx, gsm, title, spots, bases, size, strategy, source, layout}]"""
    res = json_data["result"]
    rows = []
    for uid in res["uids"]:
        r = res[uid]
        exp = r.get("expxml", "")
        runs_xml = r.get("runs", "")
        run_m = re.search(r'<Run acc="([^"]+)"[^>]*total_spots="([^"]+)"[^>]*total_bases="([^"]+)"', runs_xml)
        if not run_m:
            continue
        run, spots, bases = run_m.group(1), run_m.group(2), run_m.group(3)
        title = re.search(r"<Title>([^<]+)</Title>", exp)
        srx = re.search(r'<Experiment acc="([^"]+)"', exp)
        size = re.search(r'total_size="([^"]+)"', exp)
        strategy = re.search(r"<LIBRARY_STRATEGY>([^<]+)</LIBRARY_STRATEGY>", exp)
        source = re.search(r"<LIBRARY_SOURCE>([^<]+)</LIBRARY_SOURCE>", exp)
        layout = re.search(r"<LIBRARY_LAYOUT>\s*<([A-Z]+)/?>", exp)
        t = title.group(1) if title else ""
        gsm = t.split(":")[0].strip() if ":" in t else ""
        rows.append({
            "Run": run, "Experiment": srx.group(1) if srx else "",
            "GSM": gsm, "Title": t,
            "spots": spots, "bases": bases,
            "size_MB": str(round(int(size.group(1)) / 1e6, 1)) if size else "",
            "LibraryStrategy": strategy.group(1) if strategy else "",
            "LibrarySource": source.group(1) if source else "",
            "LibraryLayout": layout.group(1) if layout else "",
        })
    return rows

RNA_STRATEGY = ("RNA", "RNA-Seq", "miRNA-Seq", "ncRNA-Seq", "ssRNA-seq", "FL-cDNA")

def is_rnaseq(row):
    s = (row.get("LibraryStrategy") or "").upper()
    return "RNA" in s


def process_gse(gse, outdir, use_subdir=False):
    """
    处理单个 GSE，返回 (ok, n_runs, msg)

    use_subdir=False（默认）: 旧行为，文件写到 outdir/ 下，文件名含 GSE 前缀
    use_subdir=True         : 新行为，写到 outdir/{GSE}/ 下，文件名不含 GSE 前缀
    """
    try:
        # [1] SOFT -> BioProject
        text = soft_text(gse)
        bp = get_bioproject(text)
        if not bp:
            return False, 0, "无 SRA 关联(纯芯片?)"
        print(f"    BioProject = {bp}")

        # [2] BioProject -> SRA UIDs
        ids, n = esearch("sra", f"{bp}[BioProject]")
        if not ids:
            print(f"    [BioProject] 字段无结果，改用裸词查询 {bp}")
            ids, n = esearch("sra", bp)
        print(f"[2/4] SRA 中 {bp} 共有 {n} 个 run")
        if not ids:
            return False, 0, "SRA 中无 run"

        # [3] esummary 分批(每批100)
        rows = []
        for i in range(0, len(ids), 100):
            batch = ids[i:i + 100]
            raw = eutils("esummary.fcgi", {"db": "sra", "id": ",".join(batch)}, retmode="json")
            rows += parse_esummary(json.loads(raw))
        print(f"[3/4] 已解析 {len(rows)}/{n} 条")
        rows.sort(key=lambda r: r["Run"])

        # [4] 输出（根据 use_subdir 决定路径和文件名）
        if use_subdir:
            write_dir = os.path.join(outdir, gse)
            os.makedirs(write_dir, exist_ok=True)
            acc_path = os.path.join(write_dir, "SRR_Acc_List_all.txt")
            csv_path = os.path.join(write_dir, "SraRunInfo.csv")
            rna_path = os.path.join(write_dir, "SRR_Acc_List_rnaseq.txt")
        else:
            write_dir = outdir
            acc_path = os.path.join(write_dir, f"SRR_Acc_List_{gse}.txt")
            csv_path = os.path.join(write_dir, f"SraRunInfo_{gse}.csv")
            rna_path = os.path.join(write_dir, f"SRR_Acc_List_{gse}_rnaseq.txt")

        rna_rows = [r for r in rows if is_rnaseq(r)]
        with open(acc_path, "w") as f:
            f.write("\n".join(r["Run"] for r in rows) + ("\n" if rows else ""))
        with open(rna_path, "w") as f:
            f.write("\n".join(r["Run"] for r in rna_rows) + ("\n" if rna_rows else ""))
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["Run", "Experiment", "GSM", "Title", "spots", "bases", "size_MB",
                                              "LibraryStrategy", "LibrarySource", "LibraryLayout"])
            w.writeheader()
            w.writerows(rows)
        print(f"[4/4] 完成! ({len(rows)} 条, RNA-seq {len(rna_rows)} 条，全部下载不过滤)")
        print(f"      写入目录: {write_dir}")
        return True, len(rows), f"OK (总{len(rows)}/RNAseq{len(rna_rows)})"
    except Exception as e:
        return False, 0, str(e)


def write_summary(outdir):
    """
    重建 GSE_SRR_summary.csv + ALL_rnaseq_SRR.txt。
    同时支持新结构（{GSE}/SraRunInfo.csv）和旧结构（SraRunInfo_{GSE}.csv）。
    """
    from collections import Counter
    rows, all_rna = [], []

    # 新结构：{GSE}/SraRunInfo.csv
    gse_subdirs = sorted(glob.glob(os.path.join(outdir, "GSE*")))
    if gse_subdirs:
        for gse_dir in gse_subdirs:
            gse = os.path.basename(gse_dir)
            f = os.path.join(gse_dir, "SraRunInfo.csv")
            if not os.path.exists(f):
                continue
            with open(f) as fh:
                rd = list(csv.DictReader(fh))
            n_all = len(rd)
            rna = [r for r in rd if is_rnaseq(r)]
            all_rna += [r["Run"] for r in rna if r.get("Run")]
            strat = Counter((r.get("LibraryStrategy") or "?") for r in rd)
            rows.append([gse, n_all, len(rna), ";".join(f"{k}:{v}" for k, v in strat.most_common())])
    else:
        # 旧结构（向下兼容）
        for f in sorted(glob.glob(os.path.join(outdir, "SraRunInfo_*.csv"))):
            gse = os.path.basename(f).replace("SraRunInfo_", "").replace(".csv", "")
            with open(f) as fh:
                rd = list(csv.DictReader(fh))
            n_all = len(rd)
            rna = [r for r in rd if is_rnaseq(r)]
            all_rna += [r["Run"] for r in rna]
            strat = Counter((r.get("LibraryStrategy") or "?") for r in rd)
            rows.append([gse, n_all, len(rna), ";".join(f"{k}:{v}" for k, v in strat.most_common())])
    rows.sort()
    summary_path = os.path.join(outdir, "GSE_SRR_summary.csv")
    with open(summary_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["GSE_ID", "n_SRR_all", "n_SRR_RNAseq", "LibraryStrategy分布"])
        w.writerows(rows)
        w.writerow(["TOTAL", sum(r[1] for r in rows), sum(r[2] for r in rows), ""])
    with open(os.path.join(outdir, "ALL_rnaseq_SRR.txt"), "w") as fh:
        fh.write("\n".join(sorted(set(all_rna))) + ("\n" if all_rna else ""))
    print(f"[汇总] {summary_path}  ({len(rows)} 个 GSE)")
    print(f"[汇总] {os.path.join(outdir, 'ALL_rnaseq_SRR.txt')}  ({len(set(all_rna))} 个唯一 RNA-seq SRR)")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    outdir = "."
    if "-o" in sys.argv:
        outdir = sys.argv[sys.argv.index("-o") + 1]
    force = "--force" in sys.argv
    os.makedirs(outdir, exist_ok=True)

    # 批量模式
    if "--csv" in sys.argv:
        csv_file = sys.argv[sys.argv.index("--csv") + 1]
        gses = []
        with open(csv_file) as f:
            reader = csv.DictReader(f)
            cols = reader.fieldnames or []
            col = "GSE_ID" if "GSE_ID" in cols else cols[0]
            for row in reader:
                g = (row.get(col) or "").strip()
                if re.fullmatch(r"GSE\d+", g, re.I):
                    gses.append(g.upper())
        gses = list(dict.fromkeys(gses))  # 去重保序
        print(f"Pipeline-0: 输入 {csv_file}，共 {len(gses)} 个 GSE (force={force})\n")
        ok, skip, fail = 0, 0, []
        for gse in gses:
            info_path = os.path.join(outdir, f"SraRunInfo_{gse}.csv")
            if not force and os.path.exists(info_path) and os.path.getsize(info_path) > 0:
                skip += 1
                print(f"=== {gse} === 已有结果，跳过 (--force 可重拉)")
                continue
            print(f"=== {gse} ===")
            success, n, msg = process_gse(gse, outdir)
            if success:
                ok += 1
                if n == 0:
                    fail.append((gse, "有SRA但无RNA-seq数据"))
            else:
                fail.append((gse, msg))
                print(f"!! {gse} 失败: {msg}")
            print()
        print(f"===== 本次: 新拉取 {ok}，跳过 {skip}，共 {len(gses)} =====")
        if fail:
            print("无RNA-seq或失败:")
            for g, m in fail:
                print(f"  {g}: {m}")
        write_summary(outdir)
        return

    # 单 GSE 模式
    gse = sys.argv[1].strip().upper()
    success, n, msg = process_gse(gse, outdir)
    if not success:
        print(f"!! {gse} 失败: {msg}")
        sys.exit(2)
    write_summary(outdir)
    print("\n下载示例:")
    print(f"  prefetch --option-file {os.path.join(outdir, f'SRR_Acc_List_{gse}_rnaseq.txt')}")

if __name__ == "__main__":
    main()
