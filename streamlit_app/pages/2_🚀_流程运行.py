"""
pages/2_🚀_流程运行.py — 流程运行模块

Tab:
  srr_list  第一步：在 UI 输入 GSE 号 → 保存到 SRR_table.txt → 自动爬取 SRR 信息
  capacity  容量预估 + 批次计划表 + 甘特图
  run       启动 Snakemake + 实时监控状态 + 日志查看
"""

import os
import sys
import subprocess
from datetime import datetime

import streamlit as st

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from streamlit_app.core import config_loader as st_cfg
from streamlit_app.core import capacity as st_cap
from streamlit_app.core import state as st_state
from streamlit_app.core import batch_ctl as bctl
from streamlit_app.core.ui_common import render_project_selector

st.set_page_config(page_title="流程运行 — RNAseq_GEO", page_icon="🚀", layout="wide")

with st.sidebar:
    st.title("🧬 RNAseq_GEO")
    st.page_link("app.py", label="🏠 首页")

st.title("🚀 流程运行")

# 项目选择器（解决子页面直达时 session_state 为空的问题）
project = st.session_state.get("project", "")
if not project:
    project = render_project_selector()
if not project:
    st.warning("请先在「⚙️ 项目配置」中创建项目")
    st.stop()

sra_info_dir = st_cfg.get_sra_info_dir(project)
species = st_cfg.get_species(project)

# ── 持久化状态：从文件加载 / 每次渲染后保存 ──────────────────────
import json as _json_persist
_UI_STATE_FILE = os.path.join(ROOT, "result", project, "ui_run_state.json")
_UI_PERSIST_KEYS = [
    "sel_batch_file_run",   # 批次文件选择
    "run_id_val_run4",      # run_id 输入框
    "run_jobs_input",       # CPU 核心数
    "auto_r_run_val",       # 自动刷新开关
    "sel_log_tab4",         # 日志文件选择
    "tail_n_tab4",          # 日志行数
    "sel_existing_rid",     # 已有 run_id 下拉
]

# 加载（仅在 key 尚未在 session_state 中时恢复，避免覆盖用户当前操作）
if os.path.exists(_UI_STATE_FILE):
    try:
        with open(_UI_STATE_FILE) as _f:
            _saved_state = _json_persist.load(_f)
        for _k, _v in _saved_state.items():
            if _k not in st.session_state:
                st.session_state[_k] = _v
    except Exception:
        pass

def _save_ui_state():
    """将当前 session_state 中的持久化 key 写入文件"""
    try:
        os.makedirs(os.path.dirname(_UI_STATE_FILE), exist_ok=True)
        _state = {k: st.session_state[k]
                  for k in _UI_PERSIST_KEYS
                  if k in st.session_state}
        with open(_UI_STATE_FILE, "w") as _f:
            _json_persist.dump(_state, _f, ensure_ascii=False)
    except Exception:
        pass
# ────────────────────────────────────────────────────────────────

tab_srr, tab_intel, tab_cap, tab_run = st.tabs([
    "① 获取 SRR 列表",
    "② GSE 智能解读",
    "③ 容量预估与批次规划",
    "④ 启动 & 监控"
])

# ─────────── Tab ①: 获取 SRR 列表 ───────────
with tab_srr:
    st.subheader("获取 SRR 列表")

    import re as _re
    import glob as _glob

    # SRR_table.txt 路径
    srr_table_path = os.path.join(sra_info_dir, "SRR_table.txt")

    # ── 读取当前 SRR_table.txt ──
    def read_srr_table(path):
        if not os.path.exists(path):
            return []
        gses = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                m = _re.search(r"(GSE\d+)", line, _re.I)
                if m:
                    gse = m.group(1).upper()
                    if gse not in gses:
                        gses.append(gse)
        return gses

    def write_srr_table(path, gse_list):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write("# SRR_table.txt — 每行一个 GSE 号（支持 # 注释）\n")
            for gse in gse_list:
                f.write(f"{gse}\n")

    current_gses = read_srr_table(srr_table_path)

    st.caption(
        f"📄 GSE 列表文件: `{srr_table_path}` "
        f"（当前 {len(current_gses)} 个 GSE）"
    )

    col_edit, col_info = st.columns([3, 2])

    with col_edit:
        st.markdown("#### ✏️ 编辑 GSE 号列表")
        st.caption("每行一个 GSE 号（如 GSE242225），支持 # 注释行，保存后自动去重")

        # text_area 展示现有列表
        init_text = "\n".join(current_gses) if current_gses else "# 在此输入 GSE 号，每行一个\nGSE242225\n"
        new_gse_text = st.text_area(
            "GSE 号列表",
            value=init_text,
            height=300,
            key="gse_input_area",
            label_visibility="collapsed"
        )

        btn1, btn2, btn3 = st.columns(3)
        if btn1.button("💾 保存列表"):
            # 解析输入框内容
            parsed = []
            for line in new_gse_text.splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                m = _re.search(r"(GSE\d+)", line, _re.I)
                if m:
                    gse = m.group(1).upper()
                    if gse not in parsed:
                        parsed.append(gse)
            write_srr_table(srr_table_path, parsed)
            st.success(f"✅ 已保存 {len(parsed)} 个 GSE 到 `{srr_table_path}`")
            st.cache_data.clear()
            st.rerun()

    with col_info:
        st.markdown("#### 📊 已获取 SRR 信息的 GSE")
        # 新目录结构：{sra_info_dir}/{GSE}/SraRunInfo.csv
        existing_gse_dirs = sorted([
            d for d in _glob.glob(os.path.join(sra_info_dir, "GSE*"))
            if os.path.isdir(d) and os.path.exists(os.path.join(d, "SraRunInfo.csv"))
        ])
        _all_fetched = [os.path.basename(d) for d in existing_gse_dirs]
        # 只统计在当前 GSE 号列表中的（移除的 GSE 不计入）
        fetched_gses = [g for g in _all_fetched if g in current_gses] if current_gses else _all_fetched

        if current_gses:
            # 对比哪些已获取，哪些待获取（以 current_gses 为准）
            pending_fetch = [g for g in current_gses if g not in fetched_gses]
            st.success(f"✅ 已获取 {len(fetched_gses)}/{len(current_gses)} 个 GSE 的 SRA 信息")
            if pending_fetch:
                st.warning(f"⏳ 待获取: {len(pending_fetch)} 个 ({', '.join(pending_fetch[:5])}{'...' if len(pending_fetch)>5 else ''})")
        elif _all_fetched:
            st.success(f"✅ 已有 {len(_all_fetched)} 个 GSE 的 SRA 信息（GSE 列表为空）")
        else:
            st.info("尚无 SRA 信息文件，请先保存 GSE 列表并运行爬取")

    st.divider()
    st.markdown("#### 🕸️ 爬取 SRR 信息")
    st.caption(
        "点击下方按钮，脚本将从 GEO FTP 下载每个 GSE 的 SOFT 文件，提取 BioProject → SRA Run 信息。\n"
        "**幂等**：已有 SraRunInfo 的 GSE 自动跳过（--force 可强制重拉）。\n"
        "输出: `SraRunInfo_<GSE>.csv` + `SRR_Acc_List_<GSE>_rnaseq.txt` + `GSE_SRR_summary.csv`"
    )

    c_force, c_run, c_single = st.columns([1, 2, 2])
    force_fetch = c_force.checkbox("--force（强制重拉）", key="force_fetch")

    if c_run.button("▶ 爬取所有待处理 GSE"):
        if not os.path.exists(srr_table_path):
            st.error(f"SRR_table.txt 不存在，请先保存 GSE 列表")
        elif not current_gses:
            st.warning("SRR_table.txt 为空，请先输入 GSE 号")
        else:
            cmd = [
                sys.executable, "workflow/scripts/00_fetch_srr.py",
                "--table", srr_table_path,
                "--outdir", sra_info_dir,
            ]
            if force_fetch:
                cmd.append("--force")
            st.info(f"执行: `{' '.join(cmd)}`")
            with st.spinner(f"爬取 {len(current_gses)} 个 GSE 的 SRR 信息（可能需要数分钟）..."):
                proc = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
            if proc.returncode == 0:
                st.success("✅ 完成")
            else:
                st.error("❌ 爬取失败")
            st.code((proc.stdout + proc.stderr)[-3000:], language=None)
            st.cache_data.clear()

    # 单个 GSE 快速测试
    with c_single.expander("🔍 单个 GSE 测试"):
        test_gse = st.text_input("GSE 号", placeholder="GSE242225", key="test_gse")
        if st.button("爬取单个 GSE", key="fetch_single"):
            if test_gse:
                cmd = [
                    sys.executable, "workflow/scripts/00_fetch_srr.py",
                    "--gse", test_gse.strip().upper(),
                    "--outdir", sra_info_dir,
                ]
                with st.spinner(f"爬取 {test_gse}..."):
                    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
                st.code((proc.stdout + proc.stderr)[-2000:], language=None)


# ─────────── Tab ②: GSE 智能解读 ───────────
with tab_intel:
    st.subheader("GSE 智能解读")
    st.caption(
        "在运行流程前，先让 AI 阅读每个 GSE 的 GEO 信息，了解研究目的、实验设计、样本分组。\n\n"
        "**流程**：① 预拉取 GEO 信息（无需 LLM API）→ ② 输入研究背景 → ③ AI 生成解读卡"
    )

    from streamlit_app.core import geo as st_geo

    @st.cache_data(ttl=120)
    def get_sra_rows_intel(sra_dir):
        return st_geo.load_all_sra_info(sra_dir)

    sra_rows_intel = get_sra_rows_intel(sra_info_dir)
    # 只展示当前 SRR_table.txt 中的 GSE（移除的 GSE 不再展示解读卡）
    _all_fetched_gses = sorted(set(r.get("GSE", "") for r in sra_rows_intel if r.get("GSE")))
    gse_list_intel = [g for g in _all_fetched_gses if g in current_gses] if current_gses else _all_fetched_gses

    if not gse_list_intel:
        st.info("请先在「① 获取 SRR 列表」Tab 中完成 SRR 信息爬取，才能进行 GSE 解读。")
    else:
        # ── 用户背景提示输入框 ──
        st.subheader("💬 研究背景提示（可选）")
        st.caption("在此告知 AI 您的研究背景，AI 会结合 GEO 信息给出更针对性的解读。")
        col_hint1, col_hint2 = st.columns(2)
        with col_hint1:
            tumor_type = st.text_input(
                "🔬 瘤种 / 研究疾病",
                placeholder="如: GBM（胶质母细胞瘤）、肺腺癌、结直肠癌...",
                key="tumor_type_hint"
            )
            omics_type = st.selectbox(
                "📊 组学类型",
                ["RNA-seq（转录组）", "ChIP-seq（染色质免疫沉淀）",
                 "ATAC-seq（染色质可及性）", "其他"],
                key="omics_type_hint"
            )
        with col_hint2:
            extra_background = st.text_area(
                "📝 其他研究背景补充（可留空）",
                placeholder="如: 关注GBM中EGFR突变、关注IDH野生型、需要配对正常脑组织等",
                height=100,
                key="extra_background_hint"
            )

        # 构建用户提示词后缀（传递给单个 GSE 解读）
        user_hint = ""
        if tumor_type:
            user_hint += f"\n研究疾病: {tumor_type}"
        if omics_type:
            user_hint += f"\n组学类型: {omics_type}"
        if extra_background:
            user_hint += f"\n研究背景补充: {extra_background}"

        st.divider()

        # ── 整体数据集 AI 问答 ──
        st.subheader("💬 整体数据集 AI 问答")
        st.caption(
"针对您所有 GSE 数据集的整体情况，向 AI 提问。AI 会综合研究背景和 SRA 信息进行回答。（需要配置 LLM API Key）"
        )

        global_q = st.text_area(
            "输入您对整体数据集的问题",
            placeholder=(
                "例如：\n"
                "• 这批数据集中有多少个是原代组织来源？哪些是细胞系？\n"
                "• 哪些 GSE 数据集最适合做 GBM 批量转录组分析？\n"
                "• 这批数据中有没有包含配对正常脑组织的数据集？\n"
                "• 哪些数据集样本数较少，可能影响统计功效？"
            ),
            height=100,
            key="global_question"
        )

        if st.button("🤖 向 AI 提问（整体分析）", key="global_qa_btn"):
            # 构建整体上下文：汇总所有已有解读卡的关键信息
            intel_dir = os.path.join("result", project, "00_data_intel")
            summary_parts = [f"研究背景: {user_hint.strip()}" if user_hint else ""]
            summary_parts.append(f"\n共 {len(gse_list_intel)} 个 GSE 数据集:")
            for gse in gse_list_intel:
                card_path = os.path.join(intel_dir, f"{gse}.json")
                if os.path.exists(card_path):
                    try:
                        import json as _json
                        with open(card_path) as f:
                            card = _json.load(f)
                        if "error" not in card:
                            summary_parts.append(
                                f"\n{gse}: purpose={card.get('purpose','')[:80]}, "
                                f"sample_source={card.get('sample_source_type','N/A')}, "
                                f"tissue={card.get('tissue_or_cell_line','N/A')}, "
                                f"tumor_match={card.get('tumor_match_score','N/A')}, "
                                f"n_samples={card.get('n_samples','N/A')}"
                            )
                    except Exception:
                        pass
                else:
                    # 无解读卡时用 SRA 基础信息
                    gse_rows = [r for r in sra_rows_intel if r.get("GSE") == gse]
                    if gse_rows:
                        summary_parts.append(f"\n{gse}: {len(gse_rows)} 个 Run（尚无 AI 解读卡）")

            full_context = "\n".join(summary_parts)
            qa_prompt = f"{full_context}\n\n用户问题: {global_q}"

            # 调用 LLM
            qa_system = "你是生信专家，请根据提供的 GEO 数据集汇总信息，简洁准确地回答用户的问题。中文回答，300字以内。"
            with st.spinner("AI 思考中..."):
                try:
                    # 使用简单的 subprocess 调用临时脚本
                    import json as _json
                    import os as _os
                    import subprocess as _sp

                    # 读取 LLM 配置
                    llm_cfg_path = _os.path.join(ROOT, "config", "llm.yaml")
                    env_path = _os.path.join(ROOT, "config", ".env")

                    # 加载 .env
                    if _os.path.exists(env_path):
                        with open(env_path) as ef:
                            for line in ef:
                                if "=" in line and not line.startswith("#"):
                                    k, _, v = line.partition("=")
                                    _os.environ.setdefault(k.strip(), v.strip())

                    api_key = _os.environ.get("LLM_API_KEY", "")
                    if not api_key or "your-api-key" in api_key:
                        st.warning("请先在「⚙️ 项目配置 → LLM 设置」中配置 API Key")
                    else:
                        import yaml
                        with open(llm_cfg_path) as f:
                            llm_cfg = yaml.safe_load(f).get("llm", {})
                        from openai import OpenAI
                        client = OpenAI(
                            api_key=api_key,
                            base_url=llm_cfg.get("api_base", "https://api.openai.com/v1")
                        )
                        resp = client.chat.completions.create(
                            model=llm_cfg.get("model", "gpt-4o-mini"),
                            messages=[
                                {"role": "system", "content": qa_system},
                                {"role": "user",   "content": qa_prompt}
                            ],
                            temperature=float(llm_cfg.get("temperature", 0.3)),
                            timeout=float(llm_cfg.get("timeout_sec", 60)),
                            max_tokens=600,
                        )
                        answer = resp.choices[0].message.content
                        st.info(f"**AI 回答:**\n\n{answer}")
                except ImportError:
                    st.error("需要安装 openai：`pip install openai`")
                except Exception as e:
                    st.error(f"调用失败: {e}")

        st.divider()

        # ── 批量预拉取 GEO 信息（仅 SOFT，不调用 LLM）──
        st.subheader("📥 批量预拉取 GEO 元数据")
        st.caption(
            "从 GEO FTP 下载每个 GSE 的 SOFT 文件，提取：研究目的、实验设计、分组结构、基本信息。\n"
            "**不调用 LLM API**，仅爬取 GEO 公开信息。之后再进行 AI 解读。"
        )

        # 检查哪些 GSE 已有解读卡
        intel_dir = os.path.join("result", project, "00_data_intel")
        cached_gses = []
        uncached_gses = []
        for gse in gse_list_intel:
            json_path = os.path.join(intel_dir, f"{gse}.json")
            if os.path.exists(json_path) and os.path.getsize(json_path) > 0:
                cached_gses.append(gse)
            else:
                uncached_gses.append(gse)

        col_stat1, col_stat2 = st.columns(2)
        col_stat1.metric("✅ 已有解读卡", len(cached_gses))
        col_stat2.metric("⏳ 待生成", len(uncached_gses))

        c_batch_force, c_batch_run = st.columns([1, 3])
        batch_force = c_batch_force.checkbox("强制重新生成", False, key="intel_batch_force")

        if c_batch_run.button(f"▶ 批量生成所有 GSE 解读卡（{len(gse_list_intel)} 个）"):
            progress_bar = st.progress(0)
            status_text  = st.empty()
            failed_intel = []
            for idx, gse in enumerate(gse_list_intel):
                status_text.text(f"处理 {gse} ({idx+1}/{len(gse_list_intel)})...")
                cmd = [
                    sys.executable, "workflow/scripts/05_llm_reader.py",
                    "--gse",        gse,
                    "--project",    project,
                    "--sra_info",   sra_info_dir,
                    "--species",    species,
                    "--output_dir", intel_dir,
                ]
                if batch_force:
                    cmd.append("--force")
                r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT, timeout=180)
                if r.returncode != 0:
                    failed_intel.append(gse)
                progress_bar.progress((idx + 1) / len(gse_list_intel))
            status_text.text("完成！")
            if failed_intel:
                st.warning(f"以下 GSE 解读失败（可能需要配置 LLM API Key）: {failed_intel}")
            else:
                st.success(f"✅ 全部 {len(gse_list_intel)} 个 GSE 解读完成")
            st.cache_data.clear()
            st.rerun()

        st.divider()

        # ── 单个 GSE 解读展示 ──
        st.subheader("📖 查看 GSE 解读卡")
        sel_gse_intel = st.selectbox("选择 GSE", gse_list_intel, key="sel_gse_intel")

        col_card, col_gen = st.columns([3, 1])
        with col_gen:
            gen_force = st.checkbox("强制刷新", False, key="intel_single_force")
            if st.button("🤖 生成此 GSE 解读卡"):
                cmd = [
                    sys.executable, "workflow/scripts/05_llm_reader.py",
                    "--gse",        sel_gse_intel,
                    "--project",    project,
                    "--sra_info",   sra_info_dir,
                    "--species",    species,
                    "--output_dir", intel_dir,
                ]
                if gen_force:
                    cmd.append("--force")
                with st.spinner(f"调用 LLM 解读 {sel_gse_intel}..."):
                    r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT, timeout=180)
                if r.returncode == 0:
                    st.success("✅ 完成")
                else:
                    st.error(r.stderr[:500])
                st.cache_data.clear()
                st.rerun()

        with col_card:
            card = st_geo.load_llm_card(project, sel_gse_intel)
            md_card = st_geo.load_llm_card_md(project, sel_gse_intel)

            if card and "error" in card:
                st.error(f"解读失败: {card['error']}")
                if card.get("offline_mode"):
                    st.info("API 不可达。请在「⚙️ 项目配置 → LLM 设置」中配置 API Key。")
            elif md_card:
                score = card.get("reusability_score", "N/A") if card else "N/A"
                score_color = "🟢" if isinstance(score, int) and score >= 70 else \
                              "🟡" if isinstance(score, int) and score >= 40 else "🔴"
                st.info(f"**可复用性评分:** {score_color} **{score}/100**")
                if user_hint:
                    st.caption(f"💡 研究背景提示已记录（将用于 AI 解读）: {user_hint.strip()}")
                st.markdown(md_card)
            else:
                # 未有解读卡时，展示 SRA 基础信息
                gse_rows = [r for r in sra_rows_intel if r.get("GSE") == sel_gse_intel]
                if gse_rows:
                    n_samples = len(gse_rows)
                    strategies = set(r.get("LibraryStrategy", "?") for r in gse_rows)
                    st.info(
                        f"**{sel_gse_intel}** — {n_samples} 个 Run | "
                        f"策略: {', '.join(strategies)}\n\n"
                        "尚无 AI 解读卡。点击右侧「生成此 GSE 解读卡」（需配置 LLM API）。"
                    )
                else:
                    st.info(f"点击「🤖 生成此 GSE 解读卡」生成解读")


# ─────────── Tab ③: 容量预估 ───────────
with tab_cap:
    st.subheader("容量预估与批次规划")

    @st.cache_data(ttl=60)
    def load_plan(proj):
        return st_cap.load_plan(proj)

    @st.cache_data(ttl=60)
    def load_plan_df(proj):
        return st_cap.load_plan_csv(proj)

    plan = load_plan(project)

    col_btn1, col_btn2 = st.columns([1, 4])
    if col_btn1.button("🔄 重新估算"):
        cmd = [
            sys.executable, "workflow/scripts/10_capacity_planner.py",
            "--sra_info",   sra_info_dir,
            "--config",     "config/config.yaml",
            "--project",    project,
            "--output_dir", f"result/{project}/00_planning",
            "--gse_list",   srr_table_path,  # 只规划当前 GSE 列表中的 GSE
        ]
        with st.spinner("运行容量规划..."):
            r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
        if r.returncode == 0:
            st.success("✅ 完成")
            st.cache_data.clear()
            st.rerun()
        else:
            st.error(r.stderr[:500])

    if not plan:
        st.info(
            f"尚无估算结果。\n\n"
            "请先完成「获取 SRR 列表」步骤，然后点击「🔄 重新估算」。\n\n"
            "如果 SraRunInfo_*.csv 已存在，可直接点击估算。"
        )
    else:
        cal = plan.get("calibration_level", "L0")
        cal_badge = {"L0": "⚠️ L0 先验 ±50%", "L1": "🟡 L1 校准 ±20%",
                     "L2": "🟢 L2 滚动 ±10%"}.get(cal, cal)
        st.caption(f"精度: {cal_badge} | 生成: {plan.get('generated_at', '')[:16]}")

        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("总样本数", plan.get("total_samples", "N/A"))
        c2.metric("SRA 估算", f"{plan.get('total_sra_gb', 0):.0f} GB")
        c3.metric("最终产物", f"{plan.get('total_final_gb', 0):.0f} GB")
        c4.metric("预估耗时", f"{plan.get('total_hours', 0):.1f} h")
        disk_pct = plan.get("disk_usage_pct", 0)
        c5.metric("磁盘占比", f"{disk_pct:.0f}%",
                  delta="⚠️" if disk_pct > 50 else "✅",
                  delta_color="inverse" if disk_pct > 50 else "off")

        if disk_pct > 50:
            st.warning(f"⚠️ 最终产物将占磁盘剩余的 {disk_pct:.0f}%，建议扩容或分批清理！")

        # 甘特图
        gantt_path = st_cap.get_gantt_html_path(project)
        if gantt_path:
            with open(gantt_path) as f:
                gantt_html = f.read()
            st.html(gantt_html)

        # 批次计划表 → 替换为 GSE 样本选取
        plan_df = load_plan_df(project)
        if plan_df is not None:
            st.subheader("🔬 各 GSE 样本选取")
            st.caption(
                "展示每个 GSE 的 SraRunInfo.csv，默认全选。"
                "取消勾选不需要的样本后点「💾 保存选择」，"
                "保留的 SRR 将写入 `SRR_Acc_List_rnaseq.txt`（供下载和分析使用）。\n"
                "**自定义批次规划会自动使用每个 GSE 实际保留的样本数。**"
            )

            import pandas as _pd_sample
            import csv as _csv_s
            import glob as _glob_sel

            # 获取所有有 SraRunInfo.csv 的 GSE 列表（只展示当前 SRR_table.txt 中的 GSE）
            _gse_dirs_all_raw = sorted([
                os.path.basename(d)
                for d in _glob_sel.glob(os.path.join(sra_info_dir, "GSE*"))
                if os.path.isdir(d) and os.path.exists(os.path.join(d, "SraRunInfo.csv"))
            ])
            _gse_dirs_all = [g for g in _gse_dirs_all_raw if g in current_gses] if current_gses else _gse_dirs_all_raw

            if not _gse_dirs_all:
                st.info("未找到 SraRunInfo.csv，请先完成「获取 SRR 列表」步骤。")
            else:
                # 展示列映射
                _COL_SHOW = ["Run", "Title", "spots", "size_MB",
                             "LibraryStrategy", "LibraryLayout"]

                # 每个 GSE 一个 expander
                for _gse_id in _gse_dirs_all:
                    _run_info_path = os.path.join(sra_info_dir, _gse_id, "SraRunInfo.csv")
                    _rnaseq_txt    = os.path.join(sra_info_dir, _gse_id, "SRR_Acc_List_rnaseq.txt")

                    # 读取 SraRunInfo.csv
                    try:
                        _df_all = _pd_sample.read_csv(_run_info_path)
                    except Exception as _e:
                        st.warning(f"{_gse_id}: 读取 SraRunInfo.csv 失败 ({_e})")
                        continue

                    # 读取已保存的选择（SRR_Acc_List_rnaseq.txt），用于恢复勾选状态
                    _saved_srrs = set()
                    if os.path.exists(_rnaseq_txt):
                        import re as _re_srr
                        with open(_rnaseq_txt) as _rtf:
                            for _l in _rtf:
                                _m = _re_srr.search(r"(SRR\d+)", _l)
                                if _m:
                                    _saved_srrs.add(_m.group(1))

                    # 默认策略：
                    #   若 SRR_Acc_List_rnaseq.txt 已存在（由 00_gse_to_srr 生成或手动保存）
                    #     → 复现之前的选择（跨账户/跨环境均可复现）
                    #   若文件不存在或为空
                    #     → 默认只选 LibraryStrategy == RNA-Seq 的样本
                    #       （避免把 ChIP-seq/ATAC-seq 等混入 RNA-Seq 分析流程）
                    if not _saved_srrs:
                        if "LibraryStrategy" in _df_all.columns:
                            _saved_srrs = set(
                                _df_all[_df_all["LibraryStrategy"].str.upper() == "RNA-SEQ"]["Run"].tolist()
                            )
                        if not _saved_srrs:
                            # 无法按策略过滤时全选
                            _saved_srrs = set(_df_all["Run"].tolist())

                    _total_n   = len(_df_all)
                    _selected_n = sum(1 for r in _df_all["Run"] if r in _saved_srrs)

                    with st.expander(
                        f"📁 {_gse_id}  —  {_selected_n}/{_total_n} 样本已选",
                        expanded=False
                    ):
                        # 仅展示关键列，并加一列"保留"复选框
                        _show_cols = [c for c in _COL_SHOW if c in _df_all.columns]
                        _df_show = _df_all[["Run"] + [c for c in _show_cols if c != "Run"]].copy()
                        _df_show.insert(0, "保留", _df_show["Run"].isin(_saved_srrs))

                        # data_editor 支持交互式勾选
                        _edited = st.data_editor(
                            _df_show,
                            key=f"sample_sel_{_gse_id}",
                            hide_index=True,
                            width="stretch",
                            height=min(400, 35 * _total_n + 40),
                            column_config={
                                "保留": st.column_config.CheckboxColumn("保留", default=True),
                                "Run":  st.column_config.TextColumn("SRR 号", disabled=True),
                                "Title": st.column_config.TextColumn("样本描述", disabled=True),
                                "spots": st.column_config.NumberColumn("Spots", disabled=True),
                                "size_MB": st.column_config.NumberColumn("大小(MB)", disabled=True),
                                "LibraryStrategy": st.column_config.TextColumn("文库策略", disabled=True),
                                "LibraryLayout": st.column_config.TextColumn("文库类型", disabled=True),
                            },
                            disabled=[c for c in _df_show.columns if c != "保留"],
                        )

                        # 保存按钮：将勾选的 SRR 写入 SRR_Acc_List_rnaseq.txt
                        _keep_srrs = _edited[_edited["保留"] == True]["Run"].tolist()
                        _col_save, _col_info = st.columns([2, 4])
                        _col_info.caption(f"当前选中 {len(_keep_srrs)} / {_total_n} 个样本")
                        if _col_save.button(f"💾 保存 {_gse_id} 的选择", key=f"save_sel_{_gse_id}"):
                            with open(_rnaseq_txt, "w") as _wf:
                                for _srr in _keep_srrs:
                                    _wf.write(f"{_srr}\n")
                            st.success(f"✅ 已保存 {len(_keep_srrs)} 个样本到 {os.path.basename(_rnaseq_txt)}")
                            st.cache_data.clear()
                            st.rerun()

            st.divider()
            st.subheader("🗂️ 自定义批次规划")
            st.caption(
                "将所有 GSE 按实际保留样本数均衡分配到指定批次数。\n"
                "样本数来源：各 GSE 的 `SRR_Acc_List_rnaseq.txt`（已选样本数）。"
            )

            # 按批次数分组算法
            try:
                import pandas as pd
                import math
                import re as _re_n

                # 构建带「实际保留样本数」的 GSE 列表
                _gse_sample_counts = []
                for _gid in _gse_dirs_all:
                    _run_info_f = os.path.join(sra_info_dir, _gid, "SraRunInfo.csv")
                    _rnaseq_f   = os.path.join(sra_info_dir, _gid, "SRR_Acc_List_rnaseq.txt")

                    # 原始总数：从 SraRunInfo.csv 行数获取
                    _n_all = 0
                    if os.path.exists(_run_info_f):
                        try:
                            _n_all = sum(1 for _ in open(_run_info_f)) - 1  # 减去表头
                            _n_all = max(_n_all, 0)
                        except Exception:
                            pass

                    # 已选样本数：从 SRR_Acc_List_rnaseq.txt 计数
                    if os.path.exists(_rnaseq_f):
                        _cnt = sum(1 for _l in open(_rnaseq_f)
                                   if _re_n.search(r"SRR\d+", _l))
                    else:
                        # fallback：从 plan_df 取原始样本数
                        _row = plan_df[plan_df.get("gse", plan_df.iloc[:, 0]) == _gid] \
                            if "gse" in plan_df.columns else pd.DataFrame()
                        _cnt = int(_row["n_samples"].iloc[0]) if not _row.empty and "n_samples" in _row.columns else _n_all

                    if _cnt > 0:
                        _t_h = 0.0
                        if "gse" in plan_df.columns and "t_gse_h" in plan_df.columns:
                            _r = plan_df[plan_df["gse"] == _gid]
                            if not _r.empty:
                                _t_h = float(_r["t_gse_h"].iloc[0])
                        _sra_gb = 0.0
                        if "gse" in plan_df.columns and "sra_gb" in plan_df.columns:
                            _r = plan_df[plan_df["gse"] == _gid]
                            if not _r.empty:
                                _sra_gb = float(_r["sra_gb"].iloc[0])
                        _gse_sample_counts.append({
                            "gse":       _gid,
                            "n_all":     _n_all if _n_all > 0 else _cnt,
                            "n_samples": _cnt,
                            "t_gse_h":   _t_h,
                            "sra_gb":    _sra_gb,
                        })

                if not _gse_sample_counts:
                    st.info("没有找到有效的 SRR_Acc_List_rnaseq.txt，请先完成「获取 SRR 列表」步骤。")
                else:
                    df_gse = pd.DataFrame(_gse_sample_counts)
                    n_gses = len(df_gse)

                    # ── 各 GSE 样本数对照表 ──
                    st.markdown("**各 GSE 样本数对照（原始总数 vs 已选数）**")
                    _summary_rows = []
                    for _, _sr in df_gse.iterrows():
                        _n_a = int(_sr.get("n_all", _sr["n_samples"]))
                        _n_s = int(_sr["n_samples"])
                        _summary_rows.append({
                            "GSE": _sr["gse"],
                            "原始总数": _n_a,
                            "已选样本数": _n_s,
                            "过滤比例": f"{_n_s/_n_a*100:.0f}%" if _n_a > 0 else "N/A",
                        })
                    st.dataframe(
                        pd.DataFrame(_summary_rows),
                        width="stretch",
                        height=min(400, 35 * n_gses + 40),
                        hide_index=True,
                    )

                    col_n, col_sort = st.columns([2, 2])
                    n_batches = col_n.number_input(
                        "目标批次数",
                        min_value=1, max_value=n_gses,
                        value=min(10, n_gses),
                        help="将所有 GSE 均衡分配到 N 个批次，每批次独立运行"
                    )
                    sort_by = col_sort.selectbox(
                        "排序依据（影响分组方式）",
                        ["按样本数从大到小", "按预估耗时从大到小", "按 GSE 名称"],
                        index=0
                    )

                    df_sort = df_gse.copy()
                    if sort_by == "按样本数从大到小":
                        df_sort = df_sort.sort_values("n_samples", ascending=False)
                    elif sort_by == "按预估耗时从大到小":
                        df_sort = df_sort.sort_values("t_gse_h", ascending=False)
                    else:
                        df_sort = df_sort.sort_values("gse")

                    # 均衡分配：轮询法（round-robin），先大后小均匀分布
                    batches = [[] for _ in range(n_batches)]
                    batch_samples = [0] * n_batches

                    for _, row in df_sort.iterrows():
                        min_idx = batch_samples.index(min(batch_samples))
                        batches[min_idx].append(row)
                        batch_samples[min_idx] += int(row["n_samples"])

                    # 展示批次表
                    st.markdown(f"**分配结果（{n_batches} 个批次，共 {n_gses} 个 GSE）**")
                    batch_rows = []
                    for i, batch in enumerate(batches):
                        if not batch:
                            continue
                        gse_list = [r["gse"] for r in batch]
                        total_samples = sum(int(r["n_samples"]) for r in batch)
                        total_hours   = sum(float(r.get("t_gse_h", 0)) for r in batch)
                        total_sra_gb  = sum(float(r.get("sra_gb", 0)) for r in batch)
                        batch_rows.append({
                            "批次": f"Batch {i+1:02d}",
                            "GSE 数量": len(gse_list),
                            "GSE 列表": ", ".join(gse_list[:5]) + ("..." if len(gse_list) > 5 else ""),
                            "总样本数": total_samples,
                            "SRA(GB)": round(total_sra_gb, 1),
                            "预估耗时(h)": round(total_hours, 1),
                        })

                    batch_df = pd.DataFrame(batch_rows)
                    st.dataframe(batch_df, width="stretch", height=min(400, 35 * n_batches + 40))

                    # 展开查看每个批次详情
                    with st.expander("📋 各批次 GSE 详细列表"):
                        for i, batch in enumerate(batches):
                            if not batch:
                                continue
                            gse_names = [r["gse"] for r in batch]
                            st.markdown(f"**Batch {i+1:02d}**（{len(gse_names)} 个 GSE，{sum(int(r['n_samples']) for r in batch)} 个样本）")
                            st.code("\n".join(gse_names), language=None)

                    # 导出批次文件
                    st.markdown("#### 💾 生成各批次 GSE 列表文件")
                    st.caption(
                        "为每个批次生成独立的 txt 文件，保存到 `workflow/resources/{species}/`。\n"
                        "在「④ 启动 & 监控」Tab 选择要跑的批次文件后启动。\n"
                        "⚠️ **SRR_table.txt**（完整列表）不会被修改。"
                    )

                    if st.button("📂 生成全部批次文件", type="primary", key="gen_all_batches"):
                        # 先删除旧的批次文件，避免残留
                        import glob as _gb_del
                        _old_batches = _gb_del.glob(os.path.join(sra_info_dir, "batch*.txt"))
                        _removed = []
                        for _ob in _old_batches:
                            try:
                                os.remove(_ob)
                                _removed.append(os.path.basename(_ob))
                            except Exception:
                                pass
                        if _removed:
                            st.caption(f"🗑️ 已清除旧批次文件: {', '.join(_removed)}")

                        generated = []
                        for i, batch in enumerate(batches):
                            if not batch:
                                continue
                            gse_names = [r["gse"] for r in batch]
                            batch_name = f"batch{i+1:02d}_{len(gse_names)}gse"
                            batch_file = os.path.join(sra_info_dir, f"{batch_name}.txt")
                            with open(batch_file, "w") as f:
                                f.write(f"# Batch {i+1:02d} — {len(gse_names)} 个 GSE\n")
                                f.write(f"# 由批次规划自动生成\n")
                                for g in gse_names:
                                    f.write(f"{g}\n")
                            generated.append(batch_file)
                        st.success(f"✅ 已生成 {len(generated)} 个批次文件")
                        for f in generated:
                            st.caption(f"`{f}`")
                        st.cache_data.clear()

                    # 展示已有批次文件
                    import glob as _gb
                    existing_batches = sorted(_gb.glob(os.path.join(sra_info_dir, "batch*.txt")))
                    if existing_batches:
                        st.markdown("**已有批次文件：**")
                        for bf in existing_batches:
                            bname = os.path.basename(bf)
                            gse_count = sum(1 for l in open(bf) if l.strip() and not l.startswith("#"))
                            st.caption(f"• `{bname}` ({gse_count} 个 GSE)")

            except ImportError:
                st.warning("批次规划需要 pandas：`pip install pandas`")


# ─────────── Tab ④: 启动 & 监控 ───────────
with tab_run:
    st.subheader("④ 启动 & 监控 Snakemake 批处理")
    st.caption(
        "通过 `run_all.sh` 后台启动 snakemake（`gse_slots=1` 串行，一批一批跑）。\n"
        "UI 关闭后任务继续运行；重新打开 UI 点「启动/续跑」即自动续跑。"
    )

    _running = bctl.is_running()
    _pids    = bctl.find_snakemake_pids()

    # 运行状态横幅
    if _running:
        st.success(f"🟢 **批处理进行中** — PID: {_pids}")
    else:
        st.info("⚪ 批处理已停止，点「▶ 启动/续跑」开始")

    # 磁盘水位
    _total_gb, _free_gb = bctl.disk_usage_gb()
    if _total_gb > 0:
        _min_free = st_cfg.load_config().get("batch", {}).get("min_free_gb", 300)
        _pct = _free_gb / _total_gb
        st.progress(_pct, text=f"磁盘剩余: {_free_gb:.0f} GB / {_total_gb:.0f} GB ({_pct*100:.1f}%)")
        if _free_gb < _min_free:
            st.error(f"⚠️ 磁盘剩余 {_free_gb:.0f} GB < 告警阈值 {_min_free} GB！建议扩容后再继续。")

    # ── 批次文件选择 & run_id ──
    st.subheader("📂 选择批次文件")

    import glob as _gb_batch
    import re as _re_batch
    _batch_dir = st_cfg.get_sra_info_dir(project) if project else "workflow/resources/homo"
    _batch_files = sorted(_gb_batch.glob(os.path.join(_batch_dir, "batch*.txt")))
    _batch_options = ["（全量 SRR_table.txt）"] + [os.path.relpath(f, ROOT) for f in _batch_files]
    _sel_batch_file_run = st.selectbox("选择批次文件", _batch_options, key="sel_batch_file_run")

    _today = datetime.now().strftime("%Y%m%d")
    if _sel_batch_file_run != "（全量 SRR_table.txt）":
        _batch_name_stem = os.path.splitext(os.path.basename(_sel_batch_file_run))[0]
        _auto_run_id = f"{_today}_{_batch_name_stem}"
        _bf_path = os.path.join(ROOT, _sel_batch_file_run)
        _b_gses = []
        with open(_bf_path) as _bf_r:
            for _bl in _bf_r:
                _bl = _bl.strip()
                if _bl and not _bl.startswith("#"):
                    _bm = _re_batch.search(r"(GSE[0-9]+)", _bl, _re_batch.I)
                    if _bm: _b_gses.append(_bm.group(1).upper())
        if _b_gses:
            st.caption(f"包含 {len(_b_gses)} 个 GSE: {', '.join(_b_gses[:5])}{'...' if len(_b_gses)>5 else ''}")
    else:
        _auto_run_id = f"{_today}_all"

    # ── run_id 选择：优先从已有记录选，也可自由输入 ──
    # 扫描 result/{project}/ 下已有的 run_id 目录（只取含 batch 关键词的，最近的排前面）
    import glob as _gb_rid
    _existing_run_ids = []
    if project:
        _existing_run_ids = sorted(
            [os.path.basename(d) for d in _gb_rid.glob(os.path.join(ROOT, "result", project, "*"))
             if os.path.isdir(d)],
            reverse=True  # 最新的排前面
        )

    # 默认值优先级：submit_params.json > session_state > 自动生成
    import json as _json_rid
    _last_rid = None
    if not _last_rid and project:
        _param_files_rid = sorted(
            _gb_rid.glob(os.path.join(ROOT, "result", project, "*", "submit_params.json")),
            key=os.path.getmtime, reverse=True
        )
        if _param_files_rid:
            try:
                with open(_param_files_rid[0]) as _f:
                    _last_rid = _json_rid.load(_f).get("run_id")
            except Exception:
                pass

    _default_run_id = (
        st.session_state.get("run_id_val_run4")      # 1. session_state（手动输入过）
        or _last_rid                                   # 2. submit_params.json 中的上次 run_id
        or _auto_run_id                                # 3. 自动生成（今天日期）
    )
    # 只有 session_state 中尚无此 key 时才设默认值，避免覆盖用户当前输入
    if "run_id_val_run4" not in st.session_state:
        st.session_state["run_id_val_run4"] = _default_run_id

    _rid_col1, _rid_col2 = st.columns([3, 2])
    with _rid_col1:
        # 不传 value=，完全由 session_state 控制（避免 "set via Session State API" 警告）
        _run_id_val = st.text_input(
            "run_id（结果子目录名，可修改）",
            key="run_id_val_run4",
            help=(
                "结果保存到 result/{project}/{run_id}/\n\n"
                "⚠️ 续跑时必须与之前的 run_id 一致，否则会开新目录重新跑！\n"
                "可在右侧下拉框选择已有的 run_id。"
            )
        )
    with _rid_col2:
        if _existing_run_ids:
            _sel_existing = st.selectbox(
                "从已有记录选择",
                ["（新建）"] + _existing_run_ids,
                key="sel_existing_rid",
                help="选择已存在的 run_id 以续跑，而不是新建目录"
            )
            if _sel_existing != "（新建）" and _sel_existing != st.session_state.get("run_id_val_run4"):
                st.session_state["run_id_val_run4"] = _sel_existing
                st.rerun()
        else:
            st.caption("尚无已有记录")

    if project and _run_id_val:
        st.caption(f"📁 结果目录: `result/{project}/{_run_id_val}/`")

    st.divider()

    # ── 核心数控制 ──
    import os as _os_run
    _sys_cpu = _os_run.cpu_count() or 1
    _jobs_cfg = st_cfg.load_config().get("batch", {}).get("jobs", 32)
    if "run_jobs_input" not in st.session_state:
        st.session_state["run_jobs_input"] = int(_jobs_cfg)
    _run_jobs_col1, _run_jobs_col2 = st.columns([2, 3])
    with _run_jobs_col1:
        _run_jobs = st.number_input(
            "本次使用 CPU 核心数 (snakemake -j)",
            min_value=1, max_value=_sys_cpu * 2,
            key="run_jobs_input",
        )
    with _run_jobs_col2:
        st.info(f"🖥️ 服务器 CPU 总核心数: **{_sys_cpu}**（逻辑核）")

    # ── 控制按钮（启动 / 停止 / 恢复 / 强行停止）──
    _bc1, _bc2, _bc3, _bc4, _bc5 = st.columns(5)

    # ▶ 启动（未运行时可用）
    if _bc1.button("▶ 启动", type="primary", disabled=_running, key="run_start"):
        _bf = None if _sel_batch_file_run == "（全量 SRR_table.txt）" else os.path.join(ROOT, _sel_batch_file_run)
        _ri = _run_id_val.strip() if _run_id_val.strip() else None
        _ok, _msg = bctl.start(batch_file=_bf, run_id=_ri, jobs=int(_run_jobs))
        if _ok:
            # 持久化本次提交参数到 JSON 文件
            import json as _json
            _cfg_snap = st_cfg.load_config()
            _params = {
                "batch_file":        _sel_batch_file_run,
                "run_id":            _ri or _run_id_val,
                "jobs":              int(_run_jobs),
                "pipeline_threads":  _cfg_snap.get("pipeline_threads", 8),
                "pipeline_parallel": _cfg_snap.get("pipeline_parallel", 4),
                "sample_chunk_size": _cfg_snap.get("sample_chunk_size", 0),
                "submitted_at":      datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                # 可复现性参数
                "strandedness":          _cfg_snap.get("strandedness", "unstranded"),
                "alignment_rate_cutoff": _cfg_snap.get("alignment_rate_cutoff", 70.0),
                "hisat2_index":          _cfg_snap.get("hisat2_index", {}),
                "anno_base":             _cfg_snap.get("anno_base", "workflow/anno"),
                "remove_duplicates":     _cfg_snap.get("remove_duplicates", True),
            }
            st.session_state["last_submit_params"] = _params
            _param_path = os.path.join(ROOT, "result", project,
                                       _ri or _run_id_val, "submit_params.json")
            os.makedirs(os.path.dirname(_param_path), exist_ok=True)
            with open(_param_path, "w") as _pf:
                _json.dump(_params, _pf, ensure_ascii=False, indent=2)
            st.success(_msg)
        else:
            st.error(_msg)
        st.cache_data.clear()
        st.rerun()

    # ⏹ 停止（优雅 SIGTERM，运行时可用）
    if _bc2.button("⏹ 停止", disabled=not _running, key="run_stop",
                   help="发送 SIGTERM，Snakemake 会等**当前 GSE 完整跑完**后才退出（可能数小时）。\n"
                        "如需立即停止，请用「🔴 强行停止」（有数据丢失风险）。"):
        with st.spinner("发送 SIGTERM，等待当前 GSE 完成后退出（可能需要较长时间）..."):
            _ok, _msg = bctl.stop(timeout=30)
        if _ok:
            st.success(f"{_msg}\n\n⚠️ 注意：Snakemake 会等当前 GSE 跑完才真正退出，状态可能延迟更新。")
        else:
            st.warning(_msg)
        st.rerun()

    # ↺ 恢复（续跑，等同于重新启动，未运行时可用）
    if _bc3.button("↺ 恢复", disabled=_running, key="run_resume"):
        _bf = None if _sel_batch_file_run == "（全量 SRR_table.txt）" else os.path.join(ROOT, _sel_batch_file_run)
        _ri = _run_id_val.strip() if _run_id_val.strip() else None
        _ok, _msg = bctl.start(batch_file=_bf, run_id=_ri, jobs=int(_run_jobs))
        if _ok:
            import json as _json
            _cfg_snap = st_cfg.load_config()
            _params = {
                "batch_file":        _sel_batch_file_run,
                "run_id":            _ri or _run_id_val,
                "jobs":              int(_run_jobs),
                "pipeline_threads":  _cfg_snap.get("pipeline_threads", 8),
                "pipeline_parallel": _cfg_snap.get("pipeline_parallel", 4),
                "sample_chunk_size": _cfg_snap.get("sample_chunk_size", 0),
                "submitted_at":      datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                # 可复现性参数
                "strandedness":          _cfg_snap.get("strandedness", "unstranded"),
                "alignment_rate_cutoff": _cfg_snap.get("alignment_rate_cutoff", 70.0),
                "hisat2_index":          _cfg_snap.get("hisat2_index", {}),
                "anno_base":             _cfg_snap.get("anno_base", "workflow/anno"),
                "remove_duplicates":     _cfg_snap.get("remove_duplicates", True),
            }
            st.session_state["last_submit_params"] = _params
            _param_path = os.path.join(ROOT, "result", project,
                                       _ri or _run_id_val, "submit_params.json")
            os.makedirs(os.path.dirname(_param_path), exist_ok=True)
            with open(_param_path, "w") as _pf:
                _json.dump(_params, _pf, ensure_ascii=False, indent=2)
            st.success(f"🔄 续跑已启动: {_msg}")
        else:
            st.error(_msg)
        st.rerun()

    # 🔴 强行停止（SIGKILL，立即杀死，运行时可用）
    if _bc4.button("🔴 强行停止", disabled=not _running, key="run_force_stop",
                   help="立即 kill snakemake 进程（SIGKILL），强停后需点「unlock」解锁再恢复"):
        _ok, _msg = bctl.force_stop()
        if _ok:
            st.warning(f"🔴 {_msg}（提示：请点下方「unlock」按钮后再恢复）")
        else:
            st.info(_msg)
        st.rerun()

    # 🔄 刷新 + 自动刷新控制
    with _bc5:
        if st.button("🔄 刷新", key="run_manual_refresh"):
            st.rerun()
        if "auto_r_run_val" not in st.session_state:
            st.session_state["auto_r_run_val"] = False
        _auto_r = st.checkbox("15s 自动刷",
                              key="auto_r_run_val",
                              help="每 15 秒自动刷新页面状态")

    # ── Unlock（强行停止后解除目录锁，替代 Dry-run）──
    with st.expander("🔓 Unlock（强行停止后必须先解锁才能恢复）"):
        st.caption("强行停止 snakemake 后，`.snakemake/` 下会留有目录锁，需运行 unlock 才能重新启动。")
        if st.button("执行 snakemake --unlock", key="unlock_btn",
                     disabled=_running,
                     help="仅在 snakemake 未运行时可执行"):
            _ok, _msg = bctl.unlock()
            if _ok:
                st.success(_msg)
            else:
                st.error(_msg)

    st.divider()

    # ── 当前提交的 Snakemake 参数 ──
    st.subheader("📋 当前提交的 Snakemake 参数")
    import json as _json_disp
    import glob as _gb_params

    # 优先从 session_state 取，其次扫描最近的 submit_params.json
    _last = st.session_state.get("last_submit_params", {})
    if not _last and project:
        # 扫描 result/{project}/*/submit_params.json，取最近修改的
        _param_files = sorted(
            _gb_params.glob(os.path.join(ROOT, "result", project, "*", "submit_params.json")),
            key=os.path.getmtime, reverse=True
        )
        if _param_files:
            try:
                with open(_param_files[0]) as _pf:
                    _last = _json_disp.load(_pf)
                st.session_state["last_submit_params"] = _last
            except Exception:
                pass

    if _last:
        import pandas as _pd_params
        _param_label = {
            "batch_file":        "批次文件",
            "run_id":            "结果目录 (run_id)",
            "jobs":              "CPU 核心数 (-j)",
            "pipeline_threads":  "每样本线程数",
            "pipeline_parallel": "数据集内并行样本数",
            "sample_chunk_size": "滚动窗口大小",
            "submitted_at":      "提交时间",
        }
        _rows = [{"参数": _param_label.get(k, k), "值": str(v)}
                 for k, v in _last.items() if v is not None]
        st.dataframe(
            _pd_params.DataFrame(_rows),
            width="stretch",
            hide_index=True,
        )
    else:
        st.info("尚未提交过 Snakemake 任务，点「▶ 启动」或「↺ 恢复」后此处会显示参数。")

    st.divider()

    # ── 批次运行记录（扫描各 batchXX_Xgse 子目录）──
    st.subheader("📦 批次运行记录")
    st.caption("扫描 `result/{project}/` 下以批次文件命名的子目录（格式: YYYYMMDD_batchXX_Xgse）")

    import glob as _gb_prog
    _result_base = os.path.join(ROOT, "result", project) if project else None
    # 只显示含 "batch" 关键词的子目录（用户手动选批次文件启动的结果）
    _batch_run_dirs = sorted([
        os.path.basename(d.rstrip("/"))
        for d in _gb_prog.glob(os.path.join(_result_base or "", "*/"))
        if os.path.isdir(d) and "batch" in os.path.basename(d.rstrip("/")).lower()
    ]) if _result_base and os.path.isdir(_result_base) else []

    if _batch_run_dirs:
        st.caption(f"共 {len(_batch_run_dirs)} 个批次记录")
        for _rid in reversed(_batch_run_dirs):
            _rpath = os.path.join(_result_base, _rid)
            _done_list = _gb_prog.glob(os.path.join(_rpath, "03_Align_Filter", "*", "*", "dataset_finished.txt"))
            _done_cnt = len(list(_done_list))
            _has_merge = os.path.exists(os.path.join(_rpath, "04_merge_matrices", "Merge_finished.txt"))
            _icon = "✅" if _has_merge else "🔵" if _done_cnt > 0 else "⬜"
            with st.expander(f"{_icon} `{_rid}` — {_done_cnt} GSE 完成 {'| ✅矩阵已生成' if _has_merge else ''}"):
                # 读对应批次文件获取 GSE 列表
                _batch_file_name = "_".join(_rid.split("_")[1:]) + ".txt"  # 去掉日期前缀
                _bf_path = os.path.join(st_cfg.get_sra_info_dir(project) if project else "workflow/resources/homo", _batch_file_name)
                if os.path.exists(_bf_path):
                    _b_gses = []
                    import re as _re_p
                    with open(_bf_path) as _bfr:
                        for _l in _bfr:
                            _l = _l.strip()
                            if _l and not _l.startswith("#"):
                                _bm = _re_p.search(r"(GSE[0-9]+)", _l, _re_p.I)
                                if _bm: _b_gses.append(_bm.group(1).upper())
                    if _b_gses:
                        _icon_map = {"done": "✅", "running": "🔵", "pending": "⬜",
                                     "partial": "🟡", "skipped": "⏭️", "downloaded": "⬇️"}
                        for _g in _b_gses:
                            _mk = os.path.join(_rpath, "03_Align_Filter", species, _g, "dataset_finished.txt")
                            if os.path.exists(_mk):
                                _s = open(_mk).read().strip().lower()
                                _st = "done" if _s == "done" else "skipped" if "skip" in _s else "partial"
                            else:
                                _hi = os.path.join(_rpath, "03_Align_Filter", species, _g, "hisat2file")
                                _st = "running" if os.path.isdir(_hi) and os.listdir(_hi) else "pending"
                            st.text(f"  {_icon_map.get(_st, '❓')} {_g} — {_st}")
                else:
                    st.caption(f"批次文件未找到: `{_batch_file_name}`")
    else:
        st.info("尚无批次运行记录。先在「③ 容量预估」生成批次文件，选择后启动。")

    # ── 日志查看 ──
    st.divider()
    st.subheader("📜 运行日志")

    import glob as _g_run
    _log_files = sorted(
        _g_run.glob(os.path.join(ROOT, "logs", "**", "*.log"), recursive=True),
        key=os.path.getmtime, reverse=True
    )[:20]
    _log_names = [os.path.relpath(f, ROOT) for f in _log_files]

    if _log_names:
        _sel_log = st.selectbox("选择日志", _log_names, key="sel_log_tab4")
        if "tail_n_tab4" not in st.session_state:
            st.session_state["tail_n_tab4"] = 100
        _tail_n  = st.slider("显示最后 N 行", 20, 500, key="tail_n_tab4")
        if st.button("🔄 刷新日志", key="refresh_log_tab4"):
            st.rerun()
        st.code(
            bctl.tail_log(os.path.join(ROOT, _sel_log), n=_tail_n),
            language="text"
        )
    else:
        st.caption("暂无日志")

    # 自动刷新（使用 st.fragment 局部刷新，避免阻塞交互）
    # 注意：st.fragment 需要 streamlit >= 1.33；若版本不支持则回退为 autorefresh 方式
    if _auto_r:
        try:
            from streamlit_autorefresh import st_autorefresh
            st_autorefresh(interval=15000, key="auto_refresh_run")
        except ImportError:
            # fallback: 短 sleep + rerun（兼容旧版本）
            import time as _t
            _t.sleep(15)
            st.rerun()

# ── 每次渲染结束时保存 UI 状态到文件 ──
_save_ui_state()
