#!/bin/bash
# run_all.sh — 全量批次运行脚本（中断后重跑同一命令即自动续跑）
#
# 用法:
#   bash run_all.sh                          # 前台运行
#   nohup bash run_all.sh > /dev/null 2>&1 & # 后台（日志写 logs/run_all_console.log）
#
# 环境变量覆盖（Streamlit UI 通过 batch_ctl.start() 设置）:
#   SNAKEMAKE_JOBS=32          → snakemake -j
#   SNAKEMAKE_GSE_SLOTS=1      → --resources gse_slots=N
#
# 续跑语义:
#   snakemake 原生：marker 存在的批次秒跳过，只运行未完成批次。
#   中断后重跑同一命令即自动续跑，无需任何额外操作。

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# ── 读取要激活的 conda 环境名（支持环境变量覆盖） ──
CONDA_ENV="${SNAKEMAKE_CONDA_ENV:-$(python3 -c "
import yaml
try:
    cfg = yaml.safe_load(open('config/config.yaml'))
    print(cfg.get('conda_env', 'RNAseq_Pipline'))
except: print('RNAseq_Pipline')
" 2>/dev/null)}"

# ── 激活 conda 环境（非交互 shell 必须显式 source） ──
CONDA_INIT="${HOME}/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_INIT" ]; then
    source "$CONDA_INIT"
    conda activate "${CONDA_ENV}" 2>/dev/null || true
elif command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate "${CONDA_ENV}" 2>/dev/null || true
fi
echo "[$(date '+%Y-%m-%d %H:%M:%S')] 使用 conda 环境: ${CONDA_ENV}"

# 检查工具版本
if ! command -v snakemake &>/dev/null; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: 未找到 snakemake，请先激活 conda 环境"
    exit 1
fi

# 显示 prefetch 版本（用于确认环境正确）
PREFETCH_VER=$(prefetch --version 2>&1 | grep -oP '(?<=prefetch : )\S+' || echo "unknown")
echo "[$(date '+%Y-%m-%d %H:%M:%S')] prefetch 版本: ${PREFETCH_VER}"
if [[ "$PREFETCH_VER" == 2.* ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ⚠️  WARNING: prefetch 2.x 可能有 TLS 问题，建议用 RNAseq_Pipline 环境的 3.x 版本"
fi

# ── 读取 batch 配置（允许环境变量覆盖） ──
JOBS="${SNAKEMAKE_JOBS:-$(python3 -c "
import yaml, sys
try:
    cfg = yaml.safe_load(open('config/config.yaml'))
    print(cfg.get('batch', {}).get('jobs', 32))
except: print(32)
" 2>/dev/null)}"

GSE_SLOTS="${SNAKEMAKE_GSE_SLOTS:-$(python3 -c "
import yaml, sys
try:
    cfg = yaml.safe_load(open('config/config.yaml'))
    print(cfg.get('batch', {}).get('gse_slots', 1))
except: print(1)
" 2>/dev/null)}"

# ── 磁盘空间检查 ──
MIN_FREE="${MIN_FREE_GB:-300}"
FREE_GB=$(df -BG --output=avail /home | tail -1 | tr -dc '0-9' || echo 999)
if [ "${FREE_GB:-0}" -lt "$MIN_FREE" ] 2>/dev/null; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ⚠️  磁盘剩余 ${FREE_GB}G < 阈值 ${MIN_FREE}G，暂停 30min 后重试..."
    sleep 1800
    exec bash "$0"   # 重新执行自身
fi

# ── 批次文件和 run_id（支持环境变量覆盖） ──
BATCH_FILE="${SNAKEMAKE_BATCH_FILE:-}"
RUN_ID="${SNAKEMAKE_RUN_ID:-}"

# 构建 --config 追加参数
EXTRA_CONFIG=""
if [ -n "${BATCH_FILE}" ] && [ -f "${BATCH_FILE}" ]; then
    EXTRA_CONFIG="${EXTRA_CONFIG} batch_file=${BATCH_FILE}"
fi
if [ -n "${RUN_ID}" ]; then
    EXTRA_CONFIG="${EXTRA_CONFIG} run_id=${RUN_ID}"
fi

echo "======================================"
echo " RNAseq_GEO 批次自动续跑"
echo " jobs=${JOBS}  gse_slots=${GSE_SLOTS}"
[ -n "${BATCH_FILE}" ] && echo " 批次文件: ${BATCH_FILE}"
[ -n "${RUN_ID}" ] && echo " run_id: ${RUN_ID}"
echo " 时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo " 磁盘剩余: ${FREE_GB}G"
echo "======================================"

# ── 启动 snakemake ──
# 日志目录包含 run_id（按批次隔离）
LOG_DIR="logs"
if [ -n "${RUN_ID}" ]; then
    # 读取 project_name
    PNAME=$(python3 -c "
import yaml
try:
    cfg = yaml.safe_load(open('config/config.yaml'))
    print(cfg.get('projects', [{}])[0].get('project_name', 'default'))
except: print('default')
" 2>/dev/null)
    LOG_DIR="logs/${PNAME}/${RUN_ID}"
fi
mkdir -p "${LOG_DIR}"
if [ -n "${EXTRA_CONFIG}" ]; then
    exec snakemake \
        -j "${JOBS}" \
        --resources "gse_slots=${GSE_SLOTS}" \
        --rerun-incomplete \
        --keep-going \
        --configfile config/config.yaml \
        --config ${EXTRA_CONFIG} \
        2>&1 | tee -a "${LOG_DIR}/run_all_console.log"
else
    exec snakemake \
        -j "${JOBS}" \
        --resources "gse_slots=${GSE_SLOTS}" \
        --rerun-incomplete \
        --keep-going \
        --configfile config/config.yaml \
        2>&1 | tee -a "${LOG_DIR}/run_all_console.log"
fi
