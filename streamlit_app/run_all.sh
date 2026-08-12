#!/bin/bash
# streamlit_app/run_all.sh — 转发到项目根目录的 run_all.sh
# 保留此文件是因为 batch_ctl.py 引用它作为启动入口。
# 实际逻辑统一在项目根目录的 run_all.sh 中维护，此处仅做转发。

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

exec bash "${PROJECT_ROOT}/run_all.sh" "$@"
