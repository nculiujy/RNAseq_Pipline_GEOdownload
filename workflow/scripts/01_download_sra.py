
import os
import sys
import subprocess
import concurrent.futures
import shutil
import glob
import time
import re

# ============================================================
# SRA 下载助手
# 优先使用 aria2c（多连接，速度大幅提升），
# 若 aria2c 不可用或下载失败则自动回退到 prefetch。
#
# aria2c 下载原理：
#   NCBI SRA 文件可通过公开 HTTPS 路径访问：
#     https://sra-pub-run-odp.s3.amazonaws.com/sra/{SRR}/{SRR}
#   使用 aria2c -x 16 -s 16 -k 1M 开启 16 条分片并发连接，
#   充分利用本地带宽，比单连接 prefetch 快 5-10 倍。
#
# URL 模板（按优先级）：
#   1. AWS S3 公开桶 (推荐，免费，无需账号)
#   2. NCBI FTP (备用)
# ============================================================

# SRA 公开 HTTPS/FTP 下载地址模板（按优先级尝试）
SRA_URL_TEMPLATES = [
    # AWS S3 公开镜像（推荐：稳定、快速）
    "https://sra-pub-run-odp.s3.amazonaws.com/sra/{srr}/{srr}",
    # NCBI FTP（备用）
    "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos/sra-pub-run-1/{srr}/{srr}.sra",
]

def which(cmd):
    """检查命令是否存在于 PATH 中"""
    return shutil.which(cmd) is not None


def get_srr_jobs(rawdata_dir, gse_only=None):
    """
    扫描 rawdata_dir 下的 SRR 列表文件，返回 (folder, srr) 任务列表。

    支持两种文件名模式（优先级依次）：
      1. {GSE}/SRR_Acc_List_rnaseq.txt  — 新结构（00_fetch_srr.py 生成）
      2. {GSE}/SRR.txt                  — 旧结构（向下兼容）

    仅读取含 SRR 号（SRR[0-9]+）的行，自动过滤注释行和空行。
    若 gse_only 不为 None，则只处理该 GSE（per-GSE 下载用）。
    """
    print(f"[Info] 扫描目录: {rawdata_dir}" + (f"（仅 {gse_only}）" if gse_only else ""))
    jobs = []
    folders = []
    # 优先 SRR_Acc_List_all.txt（不过滤策略，下载所有 SRR）
    # 向下兼容旧文件名
    SRR_FILE_CANDIDATES = ["SRR_Acc_List_all.txt", "SRR_Acc_List_rnaseq.txt", "SRR.txt"]

    # 直接扫描一级 GSE 子目录（不递归到更深层）
    for entry in sorted(os.listdir(rawdata_dir)):
        gse_dir = os.path.join(rawdata_dir, entry)
        if not os.path.isdir(gse_dir) or not entry.upper().startswith("GSE"):
            continue
        # per-GSE 过滤
        if gse_only and entry.upper() != gse_only:
            continue

        # 找到第一个存在的 SRR 列表文件
        srr_file = None
        for fname in SRR_FILE_CANDIDATES:
            candidate = os.path.join(gse_dir, fname)
            if os.path.exists(candidate):
                srr_file = candidate
                break

        if not srr_file:
            continue

        # 提取 SRR 号
        srr_list = []
        with open(srr_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                m = re.search(r"(SRR\d+)", line)
                if m:
                    srr_list.append(m.group(1))

        if not srr_list:
            print(f"[Info] {entry}: {os.path.basename(srr_file)} 无有效 SRR 号，跳过")
            continue

        folder = entry  # 直接用 GSE 名作为子目录

        print(f"[Info] {entry}: 读取 {os.path.basename(srr_file)}，含 {len(srr_list)} 个 SRR")
        folders.append(folder)
        jobs += [(folder, srr) for srr in srr_list]

    print(f"[Info] 需下载的总 Jobs: {len(jobs)}，涉及 GSE: {folders}")
    return jobs, folders


# ------------------------------------------------------------------ #
#  aria2c 下载（多连接分片，速度比 prefetch 快 5-10 倍）              #
# ------------------------------------------------------------------ #

def _build_aria2c_cmd(url, out_path, connections=16, split=16, chunk_size="1M",
                      timeout=300, retry=5, proxy=None):
    """
    构建 aria2c 命令行列表。
    参数说明：
      connections  : -x  每个服务器的最大连接数（默认 16）
      split        : -s  将文件分成多少段并行下载（默认 16）
      chunk_size   : -k  每段最小大小（默认 1M）
      timeout      : --connect-timeout / --timeout（秒）
      retry        : --max-tries 最大重试次数
      proxy        : HTTP 代理（如 "http://127.0.0.1:7890"）
    """
    out_dir  = os.path.dirname(out_path)
    out_file = os.path.basename(out_path)
    cmd = [
        "aria2c",
        "-x", str(connections),     # 多连接并行
        "-s", str(split),           # 分片数量
        "-k", chunk_size,           # 最小分片大小
        "--connect-timeout", str(timeout),
        "--timeout", str(timeout),
        "--max-tries", str(retry),
        "--retry-wait", "5",        # 重试间隔（秒）
        "--file-allocation=none",   # 禁止预分配（加快启动）
        "--auto-file-renaming=false",
        "--allow-overwrite=true",
        "-d", out_dir,
        "-o", out_file,
    ]
    if proxy:
        cmd += ["--all-proxy", proxy]
    cmd.append(url)
    return cmd


def download_with_aria2c(srr, out_path, connections=16, split=16, chunk_size="1M",
                          timeout=300, retry=5, proxy=None):
    """
    用 aria2c 尝试各 URL 模板下载 SRR。
    成功返回 True，全部失败返回 False。
    """
    for tmpl in SRA_URL_TEMPLATES:
        url = tmpl.format(srr=srr)
        cmd = _build_aria2c_cmd(url, out_path, connections=connections,
                                split=split, chunk_size=chunk_size,
                                timeout=timeout, retry=retry, proxy=proxy)
        print(f"[aria2c] 尝试 URL: {url}")
        print(f"[aria2c] 执行命令: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, check=False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT,
                                    text=True)
            print(result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout)
            if result.returncode == 0 and os.path.exists(out_path) and os.path.getsize(out_path) > 0:
                return True
            else:
                print(f"[aria2c] URL 失败（returncode={result.returncode}），尝试下一个...")
        except Exception as e:
            print(f"[aria2c] 异常: {e}，尝试下一个...")
        # 清理残留的临时文件（.aria2 控制文件）
        for tmp in [out_path + ".aria2", out_path]:
            if os.path.exists(tmp) and (not out_path == tmp or os.path.getsize(tmp) == 0):
                try:
                    os.remove(tmp)
                except OSError:
                    pass
    return False


# ------------------------------------------------------------------ #
#  prefetch 下载（原逻辑保留，作为 fallback）                         #
# ------------------------------------------------------------------ #

def download_with_prefetch(srr, target_dir, sra_path):
    """
    使用 prefetch 下载单个 SRR，返回是否成功。
    """
    cmd = [
        "prefetch",
        "--max-size", "30GB",
        "-O", target_dir,
        srr
    ]
    print(f"[prefetch] 执行命令: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                text=True)
        print(result.stdout)

        # prefetch 可能创建 srr 子文件夹，处理路径
        srr_dir = os.path.join(target_dir, srr)
        srr_file_in_dir = os.path.join(srr_dir, f"{srr}.sra")
        if os.path.exists(srr_file_in_dir):
            shutil.move(srr_file_in_dir, sra_path)
            try:
                os.rmdir(srr_dir)
            except OSError:
                pass

        return os.path.exists(sra_path) and os.path.getsize(sra_path) > 0
    except Exception as e:
        print(f"[prefetch] 异常: {e}")
        return False


# ------------------------------------------------------------------ #
#  SRA 完整性校验（vdb-validate）                                     #
# ------------------------------------------------------------------ #

def validate_sra(sra_path, srr=""):
    """
    使用 vdb-validate 校验 .sra 文件完整性。
    返回 True（通过）或 False（校验失败或工具不可用）。
    如果 vdb-validate 不在 PATH 中，打印警告但返回 True（降级为不校验）。
    """
    if not which("vdb-validate"):
        print(f"[Validate] vdb-validate 不在 PATH 中，跳过完整性校验（建议安装 sra-tools>=3.0）")
        return True
    try:
        result = subprocess.run(
            ["vdb-validate", sra_path],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode == 0:
            print(f"[Validate] {srr}: ✅ vdb-validate 通过")
            return True
        else:
            err = (result.stdout + result.stderr).strip()[:500]
            print(f"[Validate] {srr}: ❌ vdb-validate 失败（exit={result.returncode}）: {err}")
            # 删除损坏文件，让重试可以重新下载
            if os.path.exists(sra_path):
                os.remove(sra_path)
                print(f"[Validate] {srr}: 已删除损坏的 .sra 文件")
            return False
    except subprocess.TimeoutExpired:
        print(f"[Validate] {srr}: vdb-validate 超时（120s），视为通过")
        return True
    except Exception as e:
        print(f"[Validate] {srr}: vdb-validate 异常: {e}，视为通过")
        return True


# ------------------------------------------------------------------ #
#  主下载入口（aria2c 优先，prefetch 兜底）                           #
# ------------------------------------------------------------------ #

def download_sra(result_dir, folder, srr,
                 aria2c_connections=16, aria2c_split=16, aria2c_chunk="1M",
                 aria2c_timeout=300, proxy=None):
    """
    下载单个 SRR 到 result_dir/folder/{srr}.sra。
    优先使用 aria2c（多连接），失败后回退到 prefetch。
    下载成功后自动执行 vdb-validate 完整性校验。
    """
    target_dir = os.path.join(result_dir, folder)
    os.makedirs(target_dir, exist_ok=True)

    sra_path = os.path.join(target_dir, f"{srr}.sra")
    if os.path.exists(sra_path) and os.path.getsize(sra_path) > 0:
        print(f"[Skip] {folder}: {srr}.sra 已存在，跳过下载")
        return True

    print(f"[Start] {folder}: 正在下载 {srr} ...")

    # ---------- 尝试 aria2c ----------
    if which("aria2c"):
        print(f"[Info] 使用 aria2c 多连接下载 {srr}")
        ok = download_with_aria2c(
            srr, sra_path,
            connections=aria2c_connections,
            split=aria2c_split,
            chunk_size=aria2c_chunk,
            timeout=aria2c_timeout,
            proxy=proxy,
        )
        if ok:
            # 下载成功后校验完整性
            if validate_sra(sra_path, srr):
                print(f"[Success] {folder}: {srr} aria2c 下载完成且校验通过")
                return True
            else:
                print(f"[Warn] {folder}: {srr} aria2c 下载的文件校验失败，回退到 prefetch")
        else:
            print(f"[Warn] {folder}: {srr} aria2c 全部 URL 失败，回退到 prefetch")
    else:
        print(f"[Info] aria2c 未安装，直接使用 prefetch（建议: conda install -c conda-forge aria2）")

    # ---------- 回退 prefetch ----------
    if which("prefetch"):
        ok = download_with_prefetch(srr, target_dir, sra_path)
        if ok:
            # prefetch 下载成功后校验完整性
            if validate_sra(sra_path, srr):
                print(f"[Success] {folder}: {srr} prefetch 下载完成且校验通过")
                return True
            else:
                print(f"[Error] {folder}: {srr} prefetch 下载的文件校验失败")
                return False
        print(f"[Error] {folder}: {srr} prefetch 也失败")
    else:
        print(f"[Error] aria2c 和 prefetch 均不可用，请检查环境")

    return False


def get_failed_jobs(result_dir, jobs):
    """返回未下载成功的 (folder, srr) 列表"""
    failed = []
    for folder, srr in jobs:
        sra_path = os.path.join(result_dir, folder, f"{srr}.sra")
        if not os.path.exists(sra_path) or os.path.getsize(sra_path) == 0:
            failed.append((folder, srr))
    return failed


def retry_failed_jobs(result_dir, failed_jobs, max_workers, **dl_kwargs):
    if not failed_jobs:
        print("[Info] 没有需要重试的任务")
        return []
    print(f"[Info] 准备重试 {len(failed_jobs)} 个失败的 SRR")

    def job_func(args):
        folder, srr = args
        print(f"[Retry] 再次尝试下载: {folder} - {srr}")
        return download_sra(result_dir, folder, srr, **dl_kwargs)

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        list(executor.map(job_func, failed_jobs))

    still_failed = []
    for folder, srr in failed_jobs:
        sra_path = os.path.join(result_dir, folder, f"{srr}.sra")
        if not os.path.exists(sra_path) or os.path.getsize(sra_path) == 0:
            still_failed.append((folder, srr))
    return still_failed


def main():
    print("=" * 40)
    print("[Step] 开始执行下载流程")

    import argparse
    parser = argparse.ArgumentParser(description="下载 SRA 文件（aria2c 加速 + prefetch 兜底）")
    parser.add_argument("rawdata_dir", help="rawdata 根目录（含 GSE 子目录和 SRR_Acc_List_*.txt）")
    parser.add_argument("result_dir",  help="下载结果目录（.sra 输出到 result_dir/{gse}/）")
    parser.add_argument("--gse",       default=None,  help="只下载指定 GSE（如 GSE242225），不传则下载全部")

    # aria2c 调优参数（也可通过环境变量覆盖）
    parser.add_argument("--connections", type=int, default=None,
                        help="aria2c 每服务器连接数（默认 16，可用环境变量 ARIA2_CONNECTIONS 覆盖）")
    parser.add_argument("--split",       type=int, default=None,
                        help="aria2c 分片数（默认 16，可用环境变量 ARIA2_SPLIT 覆盖）")
    parser.add_argument("--chunk",       default=None,
                        help="aria2c 最小分片大小，如 1M/2M（默认 1M，可用环境变量 ARIA2_CHUNK 覆盖）")
    parser.add_argument("--timeout",     type=int, default=None,
                        help="aria2c 连接超时秒数（默认 300）")
    parser.add_argument("--proxy",       default=None,
                        help="HTTP 代理地址（如 http://127.0.0.1:7890），也可用 ALL_PROXY 环境变量")
    args = parser.parse_args()

    rawdata_dir = args.rawdata_dir
    result_dir  = args.result_dir
    gse_only    = args.gse.strip().upper() if args.gse else None

    # 环境变量优先，命令行参数次之，再用默认值
    aria2c_connections = (args.connections
                          or int(os.environ.get("ARIA2_CONNECTIONS", "16")))
    aria2c_split       = (args.split
                          or int(os.environ.get("ARIA2_SPLIT", "16")))
    aria2c_chunk       = (args.chunk
                          or os.environ.get("ARIA2_CHUNK", "1M"))
    aria2c_timeout     = (args.timeout
                          or int(os.environ.get("ARIA2_TIMEOUT", "300")))
    proxy              = (args.proxy
                          or os.environ.get("ALL_PROXY")
                          or os.environ.get("HTTP_PROXY"))

    # 并发线程数（从环境变量读取，Snakemake 通过 resources 传入）
    max_workers = int(os.environ.get("DOWNLOAD_THREADS", "8"))

    if gse_only:
        print(f"[Info] 仅下载 GSE: {gse_only}")
    print(f"[Info] 原始数据配置目录: {rawdata_dir}")
    print(f"[Info] 结果输出目录: {result_dir}")
    print(f"[Info] aria2c 参数: connections={aria2c_connections}, split={aria2c_split}, "
          f"chunk={aria2c_chunk}, timeout={aria2c_timeout}s")
    print(f"[Info] 并发任务数: {max_workers}")
    if proxy:
        print(f"[Info] 使用代理: {proxy}")

    # aria2c 可用性检查
    if which("aria2c"):
        print("[Info] aria2c 已检测到，将优先使用多连接下载")
    else:
        print("[Warn] aria2c 未安装，将使用 prefetch（速度较慢）")
        print("[Warn] 安装建议: conda install -c conda-forge aria2")

    if not os.path.exists(rawdata_dir):
        print(f"[Error] 目录不存在: {rawdata_dir}")
        sys.exit(1)

    os.makedirs(result_dir, exist_ok=True)

    jobs, folders = get_srr_jobs(rawdata_dir, gse_only=gse_only)
    if not jobs:
        print("[Info] 没有需要下载的 SRR，退出")
        sys.exit(0)

    print("[Step] 进入下载阶段")

    dl_kwargs = dict(
        aria2c_connections=aria2c_connections,
        aria2c_split=aria2c_split,
        aria2c_chunk=aria2c_chunk,
        aria2c_timeout=aria2c_timeout,
        proxy=proxy,
    )

    def job_func(args):
        folder, srr = args
        return download_sra(result_dir, folder, srr, **dl_kwargs)

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        list(executor.map(job_func, jobs))

    # === 自动重试（最多 2 次）===
    retry_times = 2
    failed_jobs = get_failed_jobs(result_dir, jobs)
    for cycle in range(retry_times):
        if not failed_jobs:
            break
        print(f"[RetryPhase] 第 {cycle + 1} 次重试 {len(failed_jobs)} 个任务")
        still_failed = retry_failed_jobs(result_dir, failed_jobs, max_workers, **dl_kwargs)
        if not still_failed:
            print("[Info] 所有重试任务已成功")
            break
        failed_jobs = still_failed
        time.sleep(5)

    # 写入完成标志文件
    still_failed_final = get_failed_jobs(result_dir, jobs)
    finished_path = os.path.join(result_dir, "finished.txt")
    if not still_failed_final:
        with open(finished_path, "w") as f:
            f.write("done\n")
        print(f"[Done] 下载全部完成，生成标志文件: {finished_path}")
    else:
        failed_list = [f"{folder}/{srr}" for folder, srr in still_failed_final]
        with open(finished_path, "w") as f:
            f.write(f"finished_with_failures: {len(still_failed_final)} SRR 下载失败\n")
            f.write("\n".join(failed_list) + "\n")
        print(f"[Warning] {len(still_failed_final)} 个 SRR 下载失败，仍写出标志文件（失败 GSE 将被 02 跳过）")

    print(f"[Info] .sra 文件保留在 {result_dir}，由 02_dataset_pipeline.py 负责 fasterq-dump 解压")
    print("=" * 40)


if __name__ == "__main__":
    main()
