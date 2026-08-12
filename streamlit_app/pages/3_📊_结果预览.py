"""
pages/3_📊_结果预览.py — 结果预览模块

Tab:
  source    数据来源（SRR 大表 + 存储分布图）
  qc        质量控制（比对率直方图 + fastp 汇总）
  matrix    表达矩阵（8 类注释切换 + 描述统计 + PCA）

注: GSE 智能解读已迁移到「🚀 流程运行 → ② GSE 智能解读」Tab。
"""

import os
import sys
import json
import glob

import streamlit as st

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from streamlit_app.core import config_loader as st_cfg
from streamlit_app.core import geo as st_geo
from streamlit_app.core.ui_common import render_project_selector

st.set_page_config(page_title="结果预览 — RNAseq_GEO", page_icon="📊", layout="wide")

with st.sidebar:
    st.title("🧬 RNAseq_GEO")
    st.page_link("app.py", label="🏠 首页")

# 项目选择器（解决子页面直达时 session_state 为空的问题）
project = st.session_state.get("project", "")
if not project:
    project = render_project_selector()
if not project:
    st.warning("请先在「⚙️ 项目配置」中创建项目")
    st.stop()

st.title("📊 结果预览")

sra_info_dir = st_cfg.get_sra_info_dir(project)
species = st_cfg.get_species(project)


# ─────────── run_id 选择器 ───────────
def _get_run_ids(project_name, result_base="result"):
    """列出 result/{project} 下所有 run_id 子目录，按修改时间降序"""
    base = os.path.join(ROOT, result_base, project_name)
    if not os.path.isdir(base):
        return []
    dirs = [
        d for d in os.listdir(base)
        if os.path.isdir(os.path.join(base, d))
        and not d.startswith("00_")  # 排除 00_final_matrices 等项目级目录
    ]
    dirs.sort(key=lambda d: os.path.getmtime(os.path.join(base, d)), reverse=True)
    return dirs


_run_ids = _get_run_ids(project)
if not _run_ids:
    st.info("尚无运行结果（`result/{project}/` 下无 run_id 目录）")
    st.stop()

_selected_run_id = st.selectbox(
    "选择 run_id（批次）",
    _run_ids,
    index=0,
    help="默认选中最近修改的批次目录"
)
run_base = os.path.join(ROOT, "result", project, _selected_run_id)

tab_src, tab_qc, tab_matrix = st.tabs([
    "🗂️ 数据来源",
    "🩺 质量控制",
    "🧬 表达矩阵",
])


# ─────────── Tab: 数据来源 ───────────
with tab_src:
    st.subheader("🗂️ SRR 数据来源")

    @st.cache_data(ttl=120)
    def load_sra_rows(sra_dir):
        return st_geo.load_all_sra_info(sra_dir)

    sra_rows = load_sra_rows(sra_info_dir)

    if not sra_rows:
        st.info(
            f"未找到 SRA 信息（目录: `{sra_info_dir}`）\n\n"
            "请先在「🚀 流程运行 → ① 获取 SRR 列表」完成 SRR 信息爬取。"
        )
    else:
        try:
            import pandas as pd
            import plotly.express as px

            df = pd.DataFrame(sra_rows)
            st.caption(f"共 {len(df)} 条 SRR 记录，来自 `{sra_info_dir}`")

            col1, col2, col3 = st.columns(3)
            gses = sorted(df["GSE"].unique().tolist()) if "GSE" in df.columns else []
            sel_gse = col1.multiselect("GSE 过滤", gses)
            strategies = sorted(df["LibraryStrategy"].dropna().unique()) \
                if "LibraryStrategy" in df.columns else []
            sel_strat = col2.multiselect("文库策略", strategies)
            layouts = sorted(df["LibraryLayout"].dropna().unique()) \
                if "LibraryLayout" in df.columns else []
            sel_layout = col3.multiselect("文库布局", layouts)

            mask = pd.Series([True] * len(df))
            if sel_gse:
                mask &= df["GSE"].isin(sel_gse)
            if sel_strat:
                mask &= df["LibraryStrategy"].isin(sel_strat)
            if sel_layout:
                mask &= df["LibraryLayout"].isin(sel_layout)
            filtered = df[mask]

            st.caption(f"显示 {len(filtered)} / {len(df)} 条")
            st.dataframe(filtered, width="stretch", height=450)

            # 各 GSE 数据量柱状图
            if "bases" in df.columns and "GSE" in df.columns:
                st.subheader("各 GSE 数据量（bases GB）")
                try:
                    gse_gb = (
                        df.groupby("GSE")["bases"]
                        .apply(lambda x: x.astype(float).sum() / 1e9)
                        .reset_index(name="bases_GB")
                    )
                    fig = px.bar(gse_gb.sort_values("bases_GB", ascending=False),
                                 x="GSE", y="bases_GB",
                                 labels={"bases_GB": "bases (GB)"},
                                 title="各 GSE 原始数据量估算")
                    fig.update_layout(height=350)
                    st.plotly_chart(fig, width="stretch")
                except Exception as e:
                    st.caption(f"图表生成失败: {e}")
        except ImportError:
            st.warning("需要 pandas 和 plotly")


# ─────────── Tab: 质量控制 ───────────
with tab_qc:
    st.subheader("🩺 比对质量控制")

    qc_csv = os.path.join(run_base, "03_Align_Filter", "alignment_quality.csv")

    if not os.path.exists(qc_csv):
        st.info(
            "alignment_quality.csv 不存在。\n\n"
            "请先完成 Snakemake 流程的 `03_filter_alignment` 步骤。"
        )
    else:
        try:
            import pandas as pd
            import plotly.express as px

            df_qc = pd.read_csv(qc_csv)
            total = len(df_qc)
            passed = (df_qc["Passed"] == "Yes").sum() if "Passed" in df_qc.columns else 0
            failed_n = total - passed

            c1, c2, c3 = st.columns(3)
            c1.metric("总样本数", total)
            c2.metric("✅ 通过（≥70%）", passed)
            c3.metric("❌ 未通过", failed_n)

            if "Alignment_Rate" in df_qc.columns:
                fig = px.histogram(
                    df_qc, x="Alignment_Rate", nbins=30,
                    title="比对率分布（红线=70% 阈值）",
                    labels={"Alignment_Rate": "比对率 (%)"},
                    color_discrete_sequence=["#1976D2"]
                )
                fig.add_vline(x=70, line_color="red", line_dash="dash",
                              annotation_text="70% 阈值",
                              annotation_position="top right")
                st.plotly_chart(fig, width="stretch")

            # 每 GSE 通过率
            if "Sample_ID" in df_qc.columns and "Passed" in df_qc.columns:
                try:
                    df_qc["GSE"] = df_qc["Sample_ID"].str.extract(r"(GSE\d+)", expand=False)
                    gse_rate = (
                        df_qc.groupby("GSE")["Passed"]
                        .apply(lambda x: (x == "Yes").mean() * 100)
                        .reset_index(name="通过率(%)")
                    )
                    fig2 = px.bar(
                        gse_rate.sort_values("通过率(%)", ascending=True),
                        x="通过率(%)", y="GSE", orientation="h",
                        title="各 GSE 样本通过率",
                        labels={"通过率(%)": "通过率 (%)"}
                    )
                    fig2.add_vline(x=70, line_color="red", line_dash="dash")
                    fig2.update_layout(height=max(300, len(gse_rate) * 18))
                    st.plotly_chart(fig2, width="stretch")
                except Exception:
                    pass

            # 失败样本详情
            if "Passed" in df_qc.columns:
                failed_df = df_qc[df_qc["Passed"] != "Yes"]
                if not failed_df.empty:
                    with st.expander(f"❌ 未通过样本（{len(failed_df)} 个）"):
                        st.dataframe(failed_df, width="stretch")

            # fastp JSON 汇总
            fastp_files = sorted(
                glob.glob(os.path.join(run_base, "02_QC", "**", "qc_reports", "*_fastp.json"),
                          recursive=True)
            )[:20]
            if fastp_files:
                st.subheader("fastp QC 摘要（前20个样本）")
                fastp_rows = []
                for fp in fastp_files:
                    srr = os.path.basename(fp).replace("_fastp.json", "")
                    try:
                        with open(fp) as f:
                            d = json.load(f)
                        s = d.get("summary", {})
                        fastp_rows.append({
                            "SRR": srr,
                            "总 reads": s.get("before_filtering", {}).get("total_reads", ""),
                            "过滤后 reads": s.get("after_filtering", {}).get("total_reads", ""),
                            "Q20 率": s.get("after_filtering", {}).get("q20_rate", ""),
                            "GC 含量": s.get("after_filtering", {}).get("gc_content", ""),
                        })
                    except Exception:
                        pass
                if fastp_rows:
                    st.dataframe(pd.DataFrame(fastp_rows), width="stretch")
            else:
                st.caption("fastp JSON 报告未找到（流程可能尚未完成 QC 步骤）")

        except ImportError:
            st.warning("需要 pandas 和 plotly：`pip install pandas plotly`")


# ─────────── Tab: 表达矩阵 ───────────
with tab_matrix:
    st.subheader("🧬 表达量矩阵")

    matrix_files = sorted(
        glob.glob(os.path.join(run_base, "04_merge_matrices", "Matrices_*", "*_matrix.csv"))
    )

    if not matrix_files:
        st.info(
            "暂无表达矩阵。\n\n"
            "请先完成 Snakemake 流程的 `04_merge_matrices` 步骤。"
        )
    else:
        try:
            import pandas as pd
            import plotly.express as px

            matrix_names = {os.path.basename(f): f for f in matrix_files}
            sel_matrix = st.selectbox("选择注释类型", list(matrix_names.keys()))
            fpath = matrix_names[sel_matrix]

            @st.cache_data(ttl=300)
            def load_matrix(path, nrows=None):
                return pd.read_csv(path, index_col=0, nrows=nrows)

            preview_rows = st.slider("预览行数", 20, 500, 100)
            df_mat = load_matrix(fpath, nrows=preview_rows)

            n_genes, n_samples = df_mat.shape
            c1, c2, c3 = st.columns(3)
            c1.metric("基因/特征数", n_genes)
            c2.metric("样本数", n_samples)
            zero_rate = (df_mat == 0).sum().sum() / (n_genes * n_samples) * 100
            c3.metric("零值率", f"{zero_rate:.1f}%")

            st.caption(f"文件: `{fpath}`")

            tab_preview, tab_stat, tab_pca = st.tabs(["预览", "描述统计", "PCA"])

            with tab_preview:
                st.dataframe(df_mat.head(50), width="stretch")

            with tab_stat:
                st.dataframe(df_mat.describe(), width="stretch")

            with tab_pca:
                try:
                    from sklearn.decomposition import PCA
                    from sklearn.preprocessing import StandardScaler

                    mat_full = load_matrix(fpath)
                    mat_filt = mat_full[mat_full.sum(axis=1) > 0]
                    if len(mat_filt) > 500:
                        top_var = mat_filt.var(axis=1).nlargest(500).index
                        mat_filt = mat_filt.loc[top_var]

                    X = mat_filt.T.fillna(0).values
                    if X.shape[0] < 2:
                        st.warning("样本数不足，无法计算 PCA")
                    else:
                        X_scaled = StandardScaler().fit_transform(X)
                        n_comp = min(3, X_scaled.shape[0], X_scaled.shape[1])
                        pca = PCA(n_components=n_comp)
                        coords = pca.fit_transform(X_scaled)
                        pca_df = pd.DataFrame(
                            coords[:, :2],
                            columns=["PC1", "PC2"],
                            index=mat_filt.columns
                        )
                        pca_df["sample"] = pca_df.index
                        var_exp = pca.explained_variance_ratio_ * 100
                        fig = px.scatter(
                            pca_df, x="PC1", y="PC2", hover_name="sample",
                            title=f"PCA（抽样 {len(mat_filt)} 特征，PC1={var_exp[0]:.1f}%, PC2={var_exp[1]:.1f}%）",
                            labels={
                                "PC1": f"PC1 ({var_exp[0]:.1f}%)",
                                "PC2": f"PC2 ({var_exp[1]:.1f}%)"
                            }
                        )
                        st.plotly_chart(fig, width="stretch")
                except ImportError:
                    st.info("PCA 需要 scikit-learn：`pip install scikit-learn`")
                except Exception as e:
                    st.error(f"PCA 计算失败: {e}")

        except ImportError:
            st.warning("需要 pandas 和 plotly：`pip install pandas plotly`")
