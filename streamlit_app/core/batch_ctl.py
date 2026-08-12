"""
streamlit_app/core/batch_ctl.py — Streamlit 与后台 snakemake 的控制/状态桥接

核心原则:
  - UI 永远不直接跑 snakemake（会卡死）；用 Popen + start_new_session 后台启动
  - 状态判断不靠内存；用 pgrep -f snakemake 探测（UI 重启后依然准确）
  - 停止用 SIGTERM（优雅）；续跑 = 再启动一次（snakemake marker 幂等）
"""

import os
import re
import glob
import signal
import subprocess
import time
from datetime import datetime

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RUN_LOG   = os.path.join(PROJECT_ROOT, "logs", "run_all_console.log")
# run_all.sh 位于 streamlit_app/ 目录下（与 batch_ctl.py 同级目录的上一级）
_STREAMLIT_APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_SH    = os.path.join(_STREAMLIT_APP, "run_all.sh")


# ──────────────────────────────────────────────
# 进程状态检测
# ──────────────────────────────────────────────

def find_snakemake_pids():
    """探测所有 snakemake 主进程 PID（不依赖 UI 内存状态）"""
    try:
        result = subprocess.run(
            ["pgrep", "-f", "snakemake"],
            capture_output=True, text=True
        )
        pids = []
        for p in result.stdout.split():
            pid = p.strip()
            if not pid:
                continue
            # 验证：过滤掉 pgrep 自身和 streamlit
            proc_cmd = f"/proc/{pid}/cmdline"
            if os.path.exists(proc_cmd):
                with open(proc_cmd, "rb") as f:
                    cmdline = f.read().decode(errors="replace")
                if "snakemake" in cmdline and "streamlit" not in cmdline:
                    pids.append(int(pid))
        return pids
    except Exception:
        return []


def is_running():
    """是否有 snakemake 进程在运行"""
    return len(find_snakemake_pids()) > 0


# ──────────────────────────────────────────────
# 启动 / 停止
# ──────────────────────────────────────────────

def start(jobs=None, gse_slots=None, conda_env=None, batch_file=None, run_id=None):
    """
    启动/续跑批处理（通过 run_all.sh）。
    - start_new_session=True：脱离 UI 进程组，UI 关闭后任务继续跑
    - jobs/gse_slots 可覆盖 run_all.sh 中从 config 读取的值（通过环境变量）
    """
    if is_running():
        return False, "snakemake 已在运行，无需重复启动"

    if not os.path.exists(RUN_SH):
        return False, f"run_all.sh 不存在: {RUN_SH}"

    os.makedirs(os.path.dirname(RUN_LOG), exist_ok=True)
    env = os.environ.copy()
    # 环境变量覆盖（run_all.sh 中通过 ${JOBS:-32} 语法支持）
    if jobs:
        env["SNAKEMAKE_JOBS"] = str(jobs)
    if gse_slots:
        env["SNAKEMAKE_GSE_SLOTS"] = str(gse_slots)
    if conda_env:
        env["SNAKEMAKE_CONDA_ENV"] = str(conda_env)
    if batch_file:
        env["SNAKEMAKE_BATCH_FILE"] = str(batch_file)
    if run_id:
        env["SNAKEMAKE_RUN_ID"] = str(run_id)

    try:
        with open(RUN_LOG, "a") as log:
            log.write(f"\n\n{'='*60}\n")
            log.write(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] 新批次启动\n")
            log.write(f"{'='*60}\n")
            p = subprocess.Popen(
                ["bash", RUN_SH],
                cwd=PROJECT_ROOT,
                stdout=log,
                stderr=subprocess.STDOUT,
                start_new_session=True,    # 脱离 UI 进程组
                env=env
            )
        return True, f"✅ 已后台启动（PID {p.pid}），日志: `logs/run_all_console.log`"
    except Exception as e:
        return False, f"启动失败: {e}"


def stop(graceful=True, timeout=60):
    """优雅停止（SIGTERM）；超时后可改 SIGKILL"""
    pids = find_snakemake_pids()
    if not pids:
        return False, "当前没有在运行的 snakemake"

    sig = signal.SIGTERM if graceful else signal.SIGKILL
    for pid in pids:
        try:
            os.kill(pid, sig)
        except ProcessLookupError:
            pass

    # 等待退出
    deadline = time.time() + timeout
    while time.time() < deadline:
        if not find_snakemake_pids():
            return True, "✅ 已优雅停止；下次启动将自动续跑"
        time.sleep(2)

    return False, f"⚠️ {timeout}s 内未退出，请尝试强制停止"


def force_stop():
    return stop(graceful=False, timeout=10)


def unlock():
    """
    运行 snakemake --unlock，解除强制停止后留下的目录锁。
    必须在 snakemake 未运行时执行。
    """
    if is_running():
        return False, "snakemake 仍在运行，请先停止后再 unlock"
    try:
        result = subprocess.run(
            ["snakemake", "--unlock", "--configfile", "config/config.yaml"],
            capture_output=True, text=True, cwd=PROJECT_ROOT, timeout=30
        )
        out = (result.stdout + result.stderr).strip()
        if result.returncode == 0:
            return True, f"✅ unlock 成功\n{out}"
        else:
            return False, f"unlock 失败（exit {result.returncode}）:\n{out}"
    except Exception as e:
        return False, f"unlock 异常: {e}"


# ──────────────────────────────────────────────
# 批次状态扫描
# ──────────────────────────────────────────────

def _parse_gse_file(filepath):
    """读取文件中的 GSE 号列表（去重保序）"""
    if not os.path.exists(filepath):
        return []
    gses = []
    with open(filepath) as f:
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


def gse_list_from_table(rawdata_dir):
    """
    读取 GSE 号列表，优先级：
      1. batch_input.txt（自定义批次，存在时优先）
      2. SRR_table.txt（完整列表）
    与 Snakefile 的 _load_gse_list_from_table 逻辑一致。
    """
    batch_input = os.path.join(rawdata_dir, "batch_input.txt")
    if os.path.exists(batch_input):
        gses = _parse_gse_file(batch_input)
        if gses:
            return gses
    return _parse_gse_file(os.path.join(rawdata_dir, "SRR_table.txt"))


def batch_status(project, species, rawdata_dir):
    """
    扫描 markers，返回每个 GSE 的状态列表。
    status:
      pending  — 尚未开始（无 .download_done 也无 dataset_finished.txt）
      downloading — 有 .download_done（下载完）但无 dataset_finished.txt 且无 hisat2 输出
      running  — 无 dataset_finished.txt 但有 hisat2file/ 输出（正在对齐/定量）
      done     — dataset_finished.txt 内容为 "done"
      skipped  — dataset_finished.txt 内容含 "skipped"
      partial  — dataset_finished.txt 内容含 "failed" 或 "warnings"
    """
    all_gses = gse_list_from_table(rawdata_dir)
    base_01  = os.path.join(PROJECT_ROOT, "result", project,
                            "01_download_sra", species, "rawdata")
    base_03  = os.path.join(PROJECT_ROOT, "result", project,
                            "03_Align_Filter", species)
    rows = []
    for gse in all_gses:
        download_done = os.path.join(base_01, gse, ".download_done")
        marker        = os.path.join(base_03, gse, "dataset_finished.txt")
        hisat2_dir    = os.path.join(base_03, gse, "hisat2file")

        if os.path.exists(marker):
            content = open(marker).read().strip().lower()
            if "skipped" in content:
                status = "skipped"
            elif "failed" in content or "warnings" in content:
                status = "partial"
            else:
                status = "done"
        elif os.path.isdir(hisat2_dir) and os.listdir(hisat2_dir):
            status = "running"
        elif os.path.exists(download_done):
            status = "downloaded"
        else:
            status = "pending"

        rows.append({"GSE": gse, "status": status})
    return rows


def tail_log(path=None, n=100):
    """读取日志文件最后 n 行"""
    log_path = path or RUN_LOG
    if not os.path.exists(log_path):
        return "(日志不存在)"
    with open(log_path, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    return "".join(lines[-n:])


def disk_usage_gb(path="/home"):
    """返回 (total_gb, free_gb)"""
    try:
        st = os.statvfs(path)
        free_gb  = st.f_bavail * st.f_frsize / 1e9
        total_gb = st.f_blocks * st.f_frsize / 1e9
        return total_gb, free_gb
    except Exception:
        return 0, 0


def rerun_gse(project, species, gse):
    """
    重跑单个 GSE：删除 dataset_finished.txt 和 .download_done，
    然后重启批处理（snakemake 会自动检测并重跑该批）
    """
    base_01 = os.path.join(PROJECT_ROOT, "result", project,
                           "01_download_sra", species, "rawdata")
    base_03 = os.path.join(PROJECT_ROOT, "result", project,
                           "03_Align_Filter", species)

    deleted = []
    for f in [os.path.join(base_03, gse, "dataset_finished.txt"),
              os.path.join(base_01, gse, ".download_done")]:
        if os.path.exists(f):
            os.remove(f)
            deleted.append(f)

    # 重启（如果已在跑则不重复启动）
    if not is_running():
        ok, msg = start()
    else:
        ok, msg = True, "snakemake 已在运行，将自动处理该 GSE"

    return ok, f"已删除 {len(deleted)} 个标志文件 → {msg}"
