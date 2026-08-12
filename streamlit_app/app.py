"""
RNAseq_GEO — Streamlit UI 首页
================================
启动: conda activate RNAseq_Pipline && streamlit run streamlit_app/app.py

导航（左侧栏自动读取 pages/ 目录）:
  ⚙️  1_⚙️_项目配置  — config.yaml 编辑 / LLM 设置 / 工具链检测
  🚀  2_🚀_流程运行  — SRR 列表获取 / 容量预估 / 运行 & 监控
  📊  3_📊_结果预览  — 数据来源 / GSE 解读 / 质量控制 / 表达矩阵
"""

import os
import sys
import glob

import streamlit as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from streamlit_app.core import config_loader as st_cfg
from streamlit_app.core import state as st_state
from streamlit_app.core import batch_ctl as bctl

st.set_page_config(
    page_title="RNAseq_GEO",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── 全局字体缩小（约 2px）──
st.markdown("""
<style>
    html, body, [class*="css"] {
        font-size: 13px !important;
    }
    .stMetric label { font-size: 12px !important; }
    .stMetric [data-testid="stMetricValue"] { font-size: 18px !important; }
    .stCaption { font-size: 11px !important; }
</style>
""", unsafe_allow_html=True)


def _get_batch_progress(project, species):
    """
    扫描当前批次（result/{project} 下最近修改的批次目录）的 GSE 完成进度。
    返回 (done_count, total_count, batch_name) 或 (0, 0, None)。
    """
    result_base = os.path.join(ROOT, "result", project)
    if not os.path.isdir(result_base):
        return 0, 0, None

    # 找最近修改的批次子目录（含 "batch" 关键词）
    batch_dirs = sorted(
        [d for d in glob.glob(os.path.join(result_base, "*"))
         if os.path.isdir(d) and "batch" in os.path.basename(d).lower()],
        key=os.path.getmtime, reverse=True
    )
    if not batch_dirs:
        return 0, 0, None

    latest = batch_dirs[0]
    batch_name = os.path.basename(latest)

    # 从批次文件名推断对应的 txt 文件
    # 目录名格式: 20260809_batch02_2gse → batch02_2gse.txt
    parts = batch_name.split("_", 1)
    batch_file_stem = parts[1] if len(parts) > 1 else batch_name
    rawdata_dir = st_cfg.get_sra_info_dir(project)
    batch_txt = os.path.join(ROOT, rawdata_dir, f"{batch_file_stem}.txt")

    total = 0
    if os.path.exists(batch_txt):
        import re as _re
        with open(batch_txt) as f:
            for line in f:
                if _re.search(r"GSE\d+", line, _re.I):
                    total += 1

    # 统计已完成 GSE 数
    done = len(glob.glob(os.path.join(
        latest, "03_Align_Filter", species, "*", "dataset_finished.txt"
    )))

    return done, total, batch_name


# ──────────────────────────────────────────────
# 侧边栏：全局项目选择（写入 session_state 供所有页面使用）
# ──────────────────────────────────────────────
with st.sidebar:
    st.title("🧬 RNAseq_GEO")

    cfg_projects = [p.get("project_name", "") for p in st_cfg.get_projects()
                    if p.get("project_name")]
    result_projects = st_state.list_projects()
    all_projects = sorted(set(cfg_projects + result_projects))

    if not all_projects:
        st.warning("config.yaml 中无 projects 配置")
        st.session_state["project"] = ""
    else:
        idx = all_projects.index(st.session_state.get("project", all_projects[0])) \
              if st.session_state.get("project") in all_projects else 0
        project = st.selectbox("项目", all_projects, index=idx)
        st.session_state["project"] = project

    st.divider()

    # ── 批处理状态（进度条 + Snakemake 状态 + 磁盘）──
    st.subheader("⚡ 批处理状态")
    _running = bctl.is_running()
    _pids = bctl.find_snakemake_pids()

    # Snakemake 运行状态
    if _running:
        st.success(f"🟢 运行中（PID: {_pids}）")
    else:
        st.caption("⚪ Snakemake 已停止")

    # 完成通知：检测最近批次是否全部完成（含 04 merge）
    _project_check = st.session_state.get("project", "")
    if _project_check and not _running:
        _species_check = st_cfg.get_species(_project_check)
        _done_c, _total_c, _bn = _get_batch_progress(_project_check, _species_check)
        if _total_c > 0 and _done_c >= _total_c:
            # 检查 04 merge 是否也完成
            _latest_batch_dir = os.path.join(ROOT, "result", _project_check)
            _batch_dirs_c = sorted(
                [d for d in glob.glob(os.path.join(_latest_batch_dir, "*"))
                 if os.path.isdir(d) and "batch" in os.path.basename(d).lower()],
                key=os.path.getmtime, reverse=True
            )
            if _batch_dirs_c:
                _merge_marker = os.path.join(_batch_dirs_c[0], "04_merge_matrices", "Merge_finished.txt")
                if os.path.exists(_merge_marker):
                    st.success(f"✅ 批次 `{_bn}` 全部完成！（{_done_c}/{_total_c} GSE + 矩阵已合并）")

    # 当前批次进度条
    _project_sidebar = st.session_state.get("project", "")
    if _project_sidebar:
        _species_sidebar = st_cfg.get_species(_project_sidebar)
        _done, _total, _batch_name = _get_batch_progress(_project_sidebar, _species_sidebar)
        if _total > 0:
            _progress_pct = _done / _total
            st.caption(f"📦 批次: `{_batch_name}`")
            st.progress(_progress_pct, text=f"{_done}/{_total} GSE 完成 ({_progress_pct*100:.0f}%)")
        elif _batch_name:
            st.caption(f"📦 批次: `{_batch_name}` — 尚无完成记录")

    # 磁盘剩余
    _total_gb, _free_gb = bctl.disk_usage_gb()
    if _total_gb > 0:
        _min_free = st_cfg.load_config().get("batch", {}).get("min_free_gb", 300)
        _color = "🔴" if _free_gb < _min_free else "🟡" if _free_gb < _min_free * 2 else "🟢"
        _used_pct = (_total_gb - _free_gb) / _total_gb * 100
        st.caption(f"{_color} 磁盘剩余: {_free_gb:.0f} GB（已用 {_used_pct:.0f}%）")
        if _free_gb < _min_free:
            st.error(f"⚠️ 磁盘剩余 {_free_gb:.0f} GB < 告警阈值 {_min_free} GB！")

    st.divider()
    st.caption("三大模块请点击左侧导航进入")


# ──────────────────────────────────────────────
# 首页：总览卡片
# ──────────────────────────────────────────────
st.title("🧬 RNAseq_GEO")
st.markdown("**GEO 公共数据集 RNA-seq 全自动分析平台**")

project = st.session_state.get("project", "")
species = st_cfg.get_species(project) if project else "homo"

col1, col2, col3 = st.columns(3)

with col1:
    st.info("### ⚙️ 项目配置\n\n编辑 `config.yaml`，设置物种/模块开关/并发参数，配置 LLM API，检测工具链版本。")

with col2:
    st.success("### 🚀 流程运行\n\n**第一步**：获取 SRR 列表（运行 00_gse_to_srr）。\n\n**后续**：容量预估、启动 Snakemake、实时监控状态。")

with col3:
    st.warning("### 📊 结果预览\n\n查看 SRR 数据来源、GSE 智能解读卡、比对率质量控制、表达量矩阵与 PCA 可视化。")

st.divider()

# ── 项目状态摘要（5 列 metrics）──
if project:
    from streamlit_app.core import capacity as st_cap

    plan  = st_cap.load_plan(project)
    state = st_state.load_pipeline_state(project)
    proj_state = st_state.get_project_state(state, project)
    gses = proj_state.get("gses", [])

    # 当前批次样本数
    _done_batch, _total_batch, _batch_name_main = _get_batch_progress(project, species)

    # 当前项目文件夹占用量（result/{project}/）
    import subprocess as _sp
    _result_dir = os.path.join(ROOT, "result", project)
    _proj_used_gb = 0.0
    if os.path.isdir(_result_dir):
        try:
            _du = _sp.run(["du", "-sb", _result_dir],
                          capture_output=True, text=True, timeout=10)
            if _du.returncode == 0:
                _proj_used_gb = int(_du.stdout.split()[0]) / 1e9
        except Exception:
            pass

    # 磁盘整体水位（用于占比计算）
    _disk_total, _disk_free = bctl.disk_usage_gb()
    _disk_used_pct_proj = (_proj_used_gb / _disk_total * 100) if _disk_total > 0 else 0

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("当前项目", project)

    if plan:
        c2.metric("总样本数", plan.get("total_samples", "N/A"),
                  help="config 中所有 GSE 的样本总数")
    else:
        c2.metric("总样本数", "未估算",
                  help="请在「🚀 流程运行」→「③ 容量预估」生成")

    # 当前批次样本数
    if _total_batch > 0:
        c3.metric("当前批次 GSE 数",
                  f"{_total_batch} 个",
                  delta=f"✅ {_done_batch} 完成" if _done_batch > 0 else None,
                  help=f"批次: {_batch_name_main}")
    else:
        c3.metric("当前批次", "无批次",
                  help="尚未运行批次，请在「🚀 流程运行」启动")

    # 当前项目文件夹占用量（不是整个磁盘）
    c4.metric("项目目录已用",
              f"{_proj_used_gb:.1f} GB" if _proj_used_gb > 0 else "0 GB",
              help=f"result/{project}/ 目录当前占用量")

    # 项目占磁盘总量的比例
    if _disk_total > 0:
        _delta_color = "inverse" if _disk_used_pct_proj > 30 else "off"
        _delta_hint = "⚠️ 偏高" if _disk_used_pct_proj > 50 else None
        c5.metric("占磁盘比例",
                  f"{_disk_used_pct_proj:.1f}%",
                  delta=_delta_hint,
                  delta_color=_delta_color,
                  help=f"项目目录占总磁盘（{_disk_total:.0f} GB）的比例")
    else:
        c5.metric("占磁盘比例", "N/A")

    # 流程进度提示
    if gses:
        done = sum(1 for g in gses if g["status"] == "done")
        failed = sum(1 for g in gses if g["status"] == "failed")
        st.caption(f"流程进度: ✅ {done} 个 GSE 完成 | ❌ {failed} 个失败 | 共 {len(gses)} 个")
    else:
        st.caption("流程尚未开始 — 进入「🚀 流程运行」启动 Snakemake")
else:
    st.info("请在左侧侧边栏选择或配置项目")
