"""
⚙️ 项目配置 — 可视化编辑 config.yaml
参照 ChIPseq_Pipline 的设计风格：widget 读写 yaml，每个参数独立控件。

Tab:
  🔧 全局参数     — 并发/线程/质控阈值/链特异性/picard 路径
  🧬 项目管理     — 项目列表增删改（物种/GSE列表/模块开关/路径）
  📊 容量规划     — planning 参数（磁盘阈值/时间先验/甘特图开关）
  🤖 LLM 设置     — API 配置 + 密钥 + 连通性测试
  🔬 工具链检测   — 工具版本一键检测
"""

import os
import sys
import glob as _glob
import subprocess

import streamlit as st

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

CONFIG_PATH  = os.path.join(ROOT, "config", "config.yaml")
LLM_CFG_PATH = os.path.join(ROOT, "config", "llm.yaml")
ENV_PATH     = os.path.join(ROOT, "config", ".env")
ANNO_DIR     = os.path.join(ROOT, "workflow", "anno")

st.set_page_config(page_title="项目配置 — RNAseq_GEO", page_icon="⚙️", layout="wide")
st.title("⚙️ 项目配置")


# ── YAML 读写 ──────────────────────────────────────────────
def load_config():
    import yaml
    with open(CONFIG_PATH, encoding="utf-8") as f:
        return yaml.safe_load(f) or {}

def save_config(data):
    import yaml
    with open(CONFIG_PATH, "w", encoding="utf-8") as f:
        yaml.dump(data, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

def load_llm_cfg():
    import yaml
    if not os.path.exists(LLM_CFG_PATH):
        return {}
    with open(LLM_CFG_PATH, encoding="utf-8") as f:
        return yaml.safe_load(f) or {}

def save_llm_cfg(data):
    import yaml
    with open(LLM_CFG_PATH, "w", encoding="utf-8") as f:
        yaml.dump(data, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

# ── 注释目录扫描 ───────────────────────────────────────────
def get_species_dirs():
    if not os.path.isdir(ANNO_DIR):
        return []
    return sorted([d for d in os.listdir(ANNO_DIR)
                   if os.path.isdir(os.path.join(ANNO_DIR, d)) and not d.startswith(".")])

def get_hisat2_dirs(species):
    """扫描 workflow/anno/{species}/ 下所有含 hisat（大小写不限）的目录"""
    sd = os.path.join(ANNO_DIR, species)
    if not os.path.isdir(sd):
        return []
    return sorted([d for d in os.listdir(sd)
                   if os.path.isdir(os.path.join(sd, d))
                   and "hisat" in d.lower()])

def get_hisat2_prefix(species, hisat2_dir):
    """返回相对于项目根目录的 HISAT2 索引前缀，或 None。"""
    folder = os.path.join(ANNO_DIR, species, hisat2_dir)
    ht2 = _glob.glob(os.path.join(folder, "*.1.ht2"))
    if not ht2:
        return None
    # basename: GRCh38.1.ht2 → GRCh38
    prefix = os.path.basename(ht2[0]).replace(".1.ht2", "")
    return os.path.relpath(os.path.join(folder, prefix), ROOT)

def get_gtf_files(species, subdir="ncRNAanno"):
    """扫描 workflow/anno/{species}/{subdir}/ 下所有 .gtf 文件"""
    folder = os.path.join(ANNO_DIR, species, subdir)
    if not os.path.isdir(folder):
        return []
    return sorted([
        os.path.relpath(f, ROOT)
        for f in _glob.glob(os.path.join(folder, "*.gtf"))
    ])


# ═══════════════════════════════════════════════
# Tab 布局
# ═══════════════════════════════════════════════
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🔧 全局参数",
    "🧬 项目管理",
    "📊 容量规划",
    "🤖 LLM 设置",
    "🔬 工具链检测",
])


# ═══════════════════════════════════════════════
# Tab 1: 全局参数
# ═══════════════════════════════════════════════
with tab1:
    cfg = load_config()

    st.subheader("并发与线程")
    col1, col2, col3 = st.columns(3)
    with col1:
        download_threads = st.number_input(
            "下载并发数 (download_threads)", 1, 64,
            int(cfg.get("download_threads", 8)),
            help="01 prefetch 同时下载的 SRR 数量"
        )
    with col2:
        pipeline_threads = st.number_input(
            "每样本线程数 (pipeline_threads)", 1, 64,
            int(cfg.get("pipeline_threads", 8)),
            help="每个样本内 fasterq/fastp/hisat2/stringtie 使用的线程数"
        )
    with col3:
        pipeline_parallel = st.number_input(
            "数据集内并行样本数 (pipeline_parallel)", 1, 32,
            int(cfg.get("pipeline_parallel", 4)),
            help="一个 GSE 内同时处理的样本数"
        )

    col4, col5, col6 = st.columns(3)
    with col4:
        sample_chunk_size = st.number_input(
            "滚动窗口大小 (sample_chunk_size)", 0, 64,
            int(cfg.get("sample_chunk_size", 4)),
            help=(
                "每次处理 x 个样本，处理完后立即删除 sra/fastq/bam，再处理下一批，防止磁盘膨胀。\n"
                "0 = 不分块（一次性处理所有样本，原始行为）\n"
                "推荐值: 4-8（根据服务器剩余磁盘空间调整）\n"
                "建议令 sample_chunk_size ≥ pipeline_parallel，否则并行数无法跑满"
            )
        )
    with col5:
        gse_slots = st.number_input(
            "同时处理 GSE 数 (gse_slots)", 1, 8,
            int(cfg.get("batch", {}).get("gse_slots", 1)),
            help=(
                "同一时刻允许并行运行的 GSE 数量（snakemake --resources gse_slots=N）。\n"
                "1 = 严格串行（默认，磁盘最省）\n"
                "2 = 两个 GSE 同时跑，CPU 利用率更高但磁盘消耗翻倍\n"
                "建议: gse_slots × pipeline_parallel × pipeline_threads ≤ 总 CPU 核心数"
            )
        )

    st.subheader("链特异性")
    strand_opts = ["unstranded", "forward", "reverse"]
    strand_labels = {
        "unstranded": "unstranded — 不设置（非链特异性文库）",
        "forward":    "forward — FR 方向（HISAT2: FR；StringTie: --fr）",
        "reverse":    "reverse — RF 方向（dUTP 方法；HISAT2: RF；StringTie: --rf）",
    }
    cur_strand = cfg.get("strandedness", "unstranded")
    strandedness = st.selectbox(
        "strandedness",
        strand_opts,
        index=strand_opts.index(cur_strand) if cur_strand in strand_opts else 0,
        format_func=lambda x: strand_labels[x],
        help="建议用 RSeQC infer_experiment.py 确认后填写"
    )

    st.subheader("质控阈值")
    align_cutoff = st.slider(
        "比对率阈值 alignment_rate_cutoff (%)",
        0, 100, int(float(cfg.get("alignment_rate_cutoff", 70))),
        help="低于此值的样本在 03_filter_alignment 中标记为 Failed"
    )

    st.subheader("工具路径")
    picard_jar = st.text_input(
        "picard_jar",
        value=cfg.get("picard_jar", "workflow/env/picard-2.18.2/picard.jar"),
        help="Picard JAR 路径（可使用项目内路径或绝对路径）"
    )

    st.subheader("注释根目录")
    col1, col2 = st.columns(2)
    with col1:
        anno_base = st.text_input(
            "anno_base",
            value=cfg.get("anno_base", "workflow/anno"),
            help="注释文件根目录（HISAT2 索引和 GTF 文件的父目录）"
        )
    with col2:
        gtf_base = st.text_input(
            "gtf_base",
            value=cfg.get("gtf_base", "workflow/anno"),
            help="04_merge_matrices 用于识别物种/注释类型的基准路径"
        )

    # ── 并发配置可行性预检 ──
    _sys_cpu = os.cpu_count() or 1
    _total_threads_needed = int(gse_slots) * int(pipeline_parallel) * int(pipeline_threads)
    _jobs_cfg = cfg.get("batch", {}).get("jobs", 32)
    if _total_threads_needed > _sys_cpu:
        st.warning(
            f"⚠️ **并发配置偏高**: gse_slots({gse_slots}) × pipeline_parallel({pipeline_parallel}) "
            f"× pipeline_threads({pipeline_threads}) = **{_total_threads_needed}** 线程，"
            f"超过系统 CPU 核心数 **{_sys_cpu}**。\n\n"
            f"建议降低参数使乘积 ≤ {_sys_cpu}。"
        )
    if _total_threads_needed > _jobs_cfg:
        st.info(
            f"ℹ️ gse_slots × parallel × threads = {_total_threads_needed} > "
            f"batch.jobs({_jobs_cfg})，"
            f"Snakemake 调度会自动串行化部分任务（实际并行度受 -j {_jobs_cfg} 限制）。"
        )

    if st.button("💾 保存全局参数", type="primary", key="save_global"):
        fresh = load_config()
        fresh["download_threads"]       = int(download_threads)
        fresh["pipeline_threads"]       = int(pipeline_threads)
        fresh["pipeline_parallel"]      = int(pipeline_parallel)
        fresh["sample_chunk_size"]      = int(sample_chunk_size)
        fresh["strandedness"]           = strandedness
        fresh["alignment_rate_cutoff"]  = float(align_cutoff)
        fresh["picard_jar"]             = picard_jar
        fresh["anno_base"]              = anno_base
        fresh["gtf_base"]               = gtf_base
        if "batch" not in fresh:
            fresh["batch"] = {}
        fresh["batch"]["gse_slots"] = int(gse_slots)
        save_config(fresh)
        st.success("✅ 全局参数已保存")


# ═══════════════════════════════════════════════
# Tab 2: 项目管理
# ═══════════════════════════════════════════════
with tab2:
    cfg = load_config()
    projects = cfg.get("projects", [])

    st.subheader("现有项目")

    MODULE_NAMES = [
        "01_download_sra",
        "02_dataset_pipeline",
        "03_filter_alignment",
        "04_merge_matrices",
    ]
    MODULE_DESC = {
        "01_download_sra":      "prefetch 下载 .sra 文件",
        "02_dataset_pipeline":  "fasterq-dump → fastp → HISAT2 → StringTie",
        "03_filter_alignment":  "扫描比对率，生成 alignment_quality.csv",
        "04_merge_matrices":    "合并 StringTie TPM 矩阵",
    }

    for i, proj in enumerate(projects):
        pname = proj.get("project_name", f"project_{i}")
        sp    = proj.get("species", "homo")
        with st.expander(f"📁 {pname}（{sp}）", expanded=False):
            col1, col2 = st.columns(2)
            with col1:
                available_sp = get_species_dirs() or ["homo", "mouse"]
                sp_idx = available_sp.index(sp) if sp in available_sp else 0
                selected_sp = st.selectbox("物种", available_sp, index=sp_idx, key=f"sp_{i}")

                st.text_input("project_name", value=pname, key=f"pname_{i}")
                st.text_input(
                    "rawdata_dir",
                    value=proj.get("rawdata_dir", f"workflow/resources/{sp}"),
                    key=f"raw_{i}",
                    help="SRR.txt 所在根目录（01 脚本递归扫描此目录）"
                )

                # HISAT2 索引选择
                st.markdown("**HISAT2 索引**")
                ht2_dirs = get_hisat2_dirs(selected_sp)
                if ht2_dirs:
                    cur_idx = proj.get("hisat2_index", {}).get(selected_sp, "")
                    cur_ht2_dir = None
                    for d in ht2_dirs:
                        if d in (cur_idx or ""):
                            cur_ht2_dir = d
                            break
                    ht2_sel = st.selectbox(
                        "HISAT2 索引目录", ht2_dirs,
                        index=ht2_dirs.index(cur_ht2_dir) if cur_ht2_dir else 0,
                        key=f"ht2dir_{i}"
                    )
                    prefix = get_hisat2_prefix(selected_sp, ht2_sel)
                    if prefix:
                        st.caption(f"✅ 索引前缀: `{prefix}`")
                    else:
                        st.warning(f"⚠️ 未在 `{ht2_sel}/` 找到 .ht2 文件")
                        prefix = cur_idx
                else:
                    st.warning(f"⚠️ `workflow/anno/{selected_sp}/` 下未找到含 'hisat' 的索引目录")
                    prefix = st.text_input("手动输入索引前缀", value=proj.get("hisat2_index", {}).get(selected_sp, ""), key=f"ht2manual_{i}")

            with col2:
                st.markdown("**模块开关**")
                modules = proj.get("modules", {})
                new_modules = {}
                for mname in MODULE_NAMES:
                    new_modules[mname] = st.checkbox(
                        mname,
                        value=modules.get(mname, False),
                        key=f"m_{i}_{mname}",
                        help=MODULE_DESC.get(mname, ""),
                    )

            btn1, btn2 = st.columns([3, 1])
            if btn1.button(f"💾 保存 {pname}", key=f"save_p_{i}"):
                fresh = load_config()
                fp = fresh["projects"][i]
                fp["project_name"]  = st.session_state[f"pname_{i}"]
                fp["species"]       = selected_sp
                fp["rawdata_dir"]   = st.session_state[f"raw_{i}"]
                if "hisat2_index" not in fp:
                    fp["hisat2_index"] = {}
                fp["hisat2_index"][selected_sp] = prefix
                fp["modules"] = new_modules
                save_config(fresh)
                st.success(f"✅ {fp['project_name']} 已保存")
                st.rerun()
            if btn2.button("🗑️ 删除", key=f"del_p_{i}", type="secondary"):
                fresh = load_config()
                del fresh["projects"][i]
                save_config(fresh)
                st.success("已删除")
                st.rerun()

    st.markdown("---")
    st.subheader("添加新项目")

    avail_sp = get_species_dirs() or ["homo", "mouse"]
    with st.form("new_project_form"):
        col1, col2 = st.columns(2)
        with col1:
            new_sp      = st.selectbox("物种", avail_sp, key="new_sp_form")
            new_pname   = st.text_input("project_name", placeholder="如 GBM_homo")
            new_rawdata = st.text_input(
                "rawdata_dir",
                value=f"workflow/resources/homo",
                help="SRR_table.txt 所在目录（01 脚本扫描此目录下的 SRR.txt）"
            )
        with col2:
            st.markdown("**模块开关（新项目默认）**")
            new_modules = {}
            for mname in MODULE_NAMES:
                default = mname in ("01_download_sra", "02_dataset_pipeline",
                                    "03_filter_alignment", "04_merge_matrices")
                new_modules[mname] = st.checkbox(mname, value=default, key=f"nm_{mname}")

        submitted = st.form_submit_button("➕ 添加项目")
        if submitted and new_pname:
            fresh = load_config()
            if "projects" not in fresh:
                fresh["projects"] = []
            fresh["projects"].append({
                "project_name": new_pname,
                "species":       new_sp,
                "rawdata_dir":   new_rawdata,
                "modules":       new_modules,
            })
            save_config(fresh)
            st.success(f"✅ {new_pname} 已添加")
            st.rerun()


# ═══════════════════════════════════════════════
# Tab 3: 容量规划参数
# ═══════════════════════════════════════════════
with tab3:
    cfg = load_config()
    p = cfg.get("planning", {})

    st.subheader("容量规划参数（10_capacity_planner.py）")
    st.caption("这些参数影响估算精度与磁盘安全阈值")

    col1, col2 = st.columns(2)
    with col1:
        safety = st.slider(
            "磁盘安全占用上限 safety_disk_ratio",
            0.0, 1.0, float(p.get("safety_disk_ratio", 0.8)), 0.05,
            help="峰值磁盘使用量 / 磁盘剩余量 超过此比例时 UI 会告警"
        )
        io_overhead = st.slider(
            "IO/网络损耗系数 io_overhead",
            0.0, 0.5, float(p.get("io_overhead", 0.15)), 0.05,
            help="估算总耗时 × (1 + io_overhead)，补偿网络/IO 延迟"
        )
    with col2:
        alpha = st.number_input(
            "时间先验系数 alpha（秒/百万 spots）",
            0.0, 60.0, float(p.get("alpha_per_m_spots", 0.12)), 0.01,
            format="%.2f",
            help="L0 先验：单样本耗时 ≈ alpha × spots_M + beta"
        )
        beta = st.number_input(
            "固定开销 beta_fixed_min（分钟）",
            0.0, 60.0, float(p.get("beta_fixed_min", 8.0)), 0.5,
            format="%.1f",
            help="每样本固定耗时（索引加载/工具启动等）"
        )

    use_gantt = st.checkbox(
        "生成甘特图 HTML use_gantt",
        value=bool(p.get("use_gantt", True))
    )

    if st.button("💾 保存容量规划参数", type="primary", key="save_planning"):
        fresh = load_config()
        fresh["planning"] = {
            "safety_disk_ratio":  float(safety),
            "io_overhead":        float(io_overhead),
            "alpha_per_m_spots":  float(alpha),
            "beta_fixed_min":     float(beta),
            "use_gantt":          use_gantt,
        }
        save_config(fresh)
        st.success("✅ 容量规划参数已保存")


# ═══════════════════════════════════════════════
# Tab 4: LLM 设置
# ═══════════════════════════════════════════════
with tab4:
    st.subheader("LLM API 配置")
    st.caption("支持 DeepSeek / 通义千问 / 本地 vLLM 等 OpenAI 兼容接口")

    lcfg = load_llm_cfg()
    llm = lcfg.get("llm", {})

    col1, col2 = st.columns(2)
    with col1:
        api_base = st.text_input(
            "api_base",
            value=llm.get("api_base", "https://api.deepseek.com/v1"),
            help="OpenAI 兼容 API 地址"
        )
        model = st.text_input(
            "model",
            value=llm.get("model", "deepseek-chat"),
            help="模型名称"
        )
        temperature = st.slider(
            "temperature", 0.0, 1.0, float(llm.get("temperature", 0.2)), 0.05
        )

    with col2:
        timeout = st.number_input(
            "timeout_sec", 10, 300, int(llm.get("timeout_sec", 120))
        )
        max_chars = st.number_input(
            "max_context_chars", 1000, 100000, int(llm.get("max_context_chars", 30000)), 1000
        )
        proxy = st.text_input(
            "proxy（可选）",
            value=llm.get("proxy", ""),
            placeholder="http://127.0.0.1:7890",
            help="服务器外网受限时使用"
        )

    if st.button("💾 保存 LLM 配置", type="primary", key="save_llm"):
        fresh_llm = load_llm_cfg()
        fresh_llm["llm"] = {
            "provider":         "openai-compatible",
            "api_base":         api_base,
            "model":            model,
            "temperature":      float(temperature),
            "timeout_sec":      int(timeout),
            "max_context_chars": int(max_chars),
            "proxy":            proxy,
        }
        save_llm_cfg(fresh_llm)
        st.success("✅ LLM 配置已保存")

    st.divider()
    st.subheader("API Key 管理")

    # 读取现有 .env（仅显示 key 名，不显示值）
    env_keys = {}
    if os.path.exists(ENV_PATH):
        with open(ENV_PATH) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    k, _, _ = line.partition("=")
                    env_keys[k.strip()] = True
    if env_keys:
        st.caption(f"config/.env 中已有 key: {', '.join(env_keys.keys())}")
    else:
        st.caption("config/.env 中尚无 key")

    new_key = st.text_input("输入 LLM_API_KEY", type="password", placeholder="sk-...")
    if st.button("💾 保存 API Key") and new_key:
        with open(ENV_PATH, "w") as f:
            f.write(f"LLM_API_KEY={new_key}\n")
        st.success("✅ API Key 已保存到 config/.env")

    st.divider()
    st.subheader("🔌 连通性测试")
    if st.button("测试 LLM 连接"):
        test_script = f"""
import os, sys
sys.path.insert(0, {repr(ROOT)})
# 加载 .env
env_path = {repr(ENV_PATH)}
if os.path.exists(env_path):
    with open(env_path) as f:
        for line in f:
            if '=' in line and not line.startswith('#'):
                k, _, v = line.partition('=')
                os.environ.setdefault(k.strip(), v.strip())
api_key = os.environ.get('LLM_API_KEY', '')
if not api_key or 'your-api-key' in api_key:
    print('ERROR: 请先配置真实 API Key')
    sys.exit(1)
try:
    from openai import OpenAI
    client = OpenAI(api_key=api_key, base_url={repr(api_base)})
    r = client.chat.completions.create(
        model={repr(model)},
        messages=[{{'role':'user','content':'Reply OK'}}],
        max_tokens=5, timeout=30)
    print('✅ 连接成功:', r.choices[0].message.content)
except Exception as e:
    print('❌ 失败:', e)
"""
        with st.spinner("测试中..."):
            r = subprocess.run([sys.executable, "-c", test_script],
                               capture_output=True, text=True, cwd=ROOT, timeout=40)
        out = (r.stdout + r.stderr).strip()
        if "✅" in out:
            st.success(out)
        else:
            st.error(out)


# ═══════════════════════════════════════════════
# Tab 5: 工具链检测
# ═══════════════════════════════════════════════
with tab5:
    st.subheader("工具链版本检测")

    # 从 config 读取 conda 环境名
    _conda_env = load_config().get("conda_env", "RNAseq_Pipline")
    st.caption(
        f"通过 `conda run -n {_conda_env}` 执行检测，"
        f"确保工具在目标环境中可用（即使 streamlit 从其他环境启动）。"
    )

    tools_info = [
        ("hisat2",        "hisat2 --version"),
        ("fastp",         "fastp --version"),
        ("stringtie",     "stringtie --version"),
        ("samtools",      "samtools --version"),
        ("fasterq-dump",  "fasterq-dump --version"),
        ("prefetch",      "prefetch --version"),
        ("snakemake",     "snakemake --version"),
        ("python",        "python --version"),
        ("java",          "java -version"),
    ]

    col_name, col_status = st.columns([2, 8])
    col_name.markdown("**工具**")
    col_status.markdown("**版本信息**")

    for name, cmd_str in tools_info:
        c1, c2 = st.columns([2, 8])
        c1.write(f"`{name}`")
        try:
            # 使用 conda run 在目标环境中执行，避免 PATH 不包含工具的误报
            full_cmd = ["conda", "run", "-n", _conda_env, "--no-banner", "bash", "-c", cmd_str]
            r = subprocess.run(full_cmd, capture_output=True, text=True, timeout=15)
            output = (r.stdout + r.stderr).strip()
            ver = output.split("\n")[0][:100] if output else ""
            if r.returncode == 0 and ver:
                c2.success(ver)
            elif ver:
                # 部分工具 --version 返回非零但输出版本号（如 java -version）
                c2.success(ver)
            else:
                c2.error(f"❌ 未找到或执行失败（exit={r.returncode}）")
        except FileNotFoundError:
            # conda 本身不在 PATH
            c2.error("❌ conda 命令不可用（请确认 conda 已安装并在 PATH 中）")
        except subprocess.TimeoutExpired:
            c2.warning("⚠️ 超时")
        except Exception as e:
            c2.warning(str(e)[:80])
