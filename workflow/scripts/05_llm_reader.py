#!/usr/bin/env python3
"""
05_llm_reader.py — GEO 数据集智能解读（LLM API）

功能:
  1. 读取 GEO SOFT 元数据 + SraRunInfo_<GSE>.csv + 可选 QC 快照
  2. 组装上下文（防超长截断）
  3. 调用 OpenAI 兼容 API 生成结构化数据解读卡
  4. 输出 <GSE>.md（人读）+ <GSE>.json（UI 渲染用）

配置:
  config/llm.yaml   → API 地址 / 模型 / 超时等
  config/.env       → LLM_API_KEY=sk-xxx

用法:
  python workflow/scripts/05_llm_reader.py \\
      --gse        GSE242225 \\
      --project    GBM_homo \\
      --sra_info   workflow/resources/homo/ \\
      --output_dir result/GBM_homo/00_data_intel \\
      [--qc_dir    result/GBM_homo/03_Align_Filter] \\
      [--force]    # 强制覆盖已有缓存

依赖（补装）:
  pip install openai python-dotenv
  或: conda install -c conda-forge openai python-dotenv
"""

import os
import sys
import csv
import json
import argparse
import re
import time
import gzip
import urllib.request
from datetime import datetime

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

try:
    from dotenv import load_dotenv
    HAS_DOTENV = True
except ImportError:
    HAS_DOTENV = False

try:
    from openai import OpenAI
    HAS_OPENAI = True
except ImportError:
    HAS_OPENAI = False


# ──────────────────────────────────────────────
# 工具函数
# ──────────────────────────────────────────────

def load_llm_config(llm_yaml_path="config/llm.yaml"):
    if not os.path.exists(llm_yaml_path):
        return {}
    with open(llm_yaml_path) as f:
        if HAS_YAML:
            cfg = yaml.safe_load(f) or {}
        else:
            cfg = {}
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and ":" in line:
                    k, _, v = line.partition(":")
                    cfg[k.strip()] = v.strip()
    return cfg.get("llm", cfg) if "llm" in cfg else cfg


def load_env(env_path="config/.env"):
    if HAS_DOTENV and os.path.exists(env_path):
        load_dotenv(env_path)
    # fallback: 手动解析
    elif os.path.exists(env_path):
        with open(env_path) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    k, _, v = line.partition("=")
                    os.environ.setdefault(k.strip(), v.strip())


def load_sra_info(sra_info_dir, gse):
    """
    读取 SRA 信息，支持新目录结构（{gse}/SraRunInfo.csv）和旧结构（SraRunInfo_{gse}.csv）
    """
    # 新结构优先
    new_path = os.path.join(sra_info_dir, gse, "SraRunInfo.csv")
    if os.path.exists(new_path):
        with open(new_path, newline="") as f:
            return list(csv.DictReader(f))
    # 旧结构兼容
    old_path = os.path.join(sra_info_dir, f"SraRunInfo_{gse}.csv")
    if os.path.exists(old_path):
        with open(old_path, newline="") as f:
            return list(csv.DictReader(f))
    return []


def load_qc_summary(qc_dir, gse, species):
    """
    扫描 qc_dir/{species}/{gse}/hisat2file/*/QC_results.log，
    提取比对率分布（最多 50 条）
    """
    pattern = os.path.join(qc_dir, species, gse, "hisat2file", "*", "QC_results.log")
    import glob
    logs = glob.glob(pattern)
    if not logs:
        return None
    lines = []
    for lpath in sorted(logs)[:50]:
        srr = os.path.basename(os.path.dirname(lpath))
        try:
            with open(lpath) as f:
                content = f.read()
            matches = re.findall(r"(\d+\.\d+)%", content)
            rate = matches[-1] if matches else "N/A"
            lines.append(f"  {srr}: {rate}%")
        except Exception:
            lines.append(f"  {srr}: 读取失败")
    return "比对率快照（前50样本）:\n" + "\n".join(lines)


def fetch_geo_soft_summary(gse):
    """
    从 GEO FTP 下载 SOFT 文件，提取 summary / overall_design / sample 表
    """
    m = re.fullmatch(r"GSE(\d+)", gse.upper())
    if not m:
        return None
    digits = m.group(1)
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
           f"GSE{digits[:3]}nnn/{gse.upper()}/soft/{gse.upper()}_family.soft.gz")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "05_llm_reader/1.0"})
        with urllib.request.urlopen(req, timeout=30) as r:
            raw = r.read()
        if raw[:2] == b"\x1f\x8b":
            text = gzip.decompress(raw).decode("utf-8", "replace")
        else:
            text = raw.decode("utf-8", "replace")
        return text
    except Exception as e:
        print(f"[WARN] SOFT 下载失败: {e}")
        return None


def parse_soft(text, max_sample_rows=50):
    """从 SOFT 文本提取 title / summary / overall_design / sample 表"""
    if not text:
        return {}
    result = {}
    for field in ["Series_title", "Series_summary", "Series_overall_design"]:
        m = re.search(rf"!{field}\s*=\s*(.+)", text)
        if m:
            result[field.replace("Series_", "")] = m.group(1).strip()
    # 样本表
    sample_blocks = re.findall(r"\^SAMPLE.*?(?=\^SAMPLE|\Z)", text, re.DOTALL)
    rows = []
    for block in sample_blocks[:max_sample_rows]:
        title_m = re.search(r"!Sample_title\s*=\s*(.+)", block)
        gsm_m   = re.search(r"!Sample_geo_accession\s*=\s*(.+)", block)
        char_m  = re.findall(r"!Sample_characteristics_ch1\s*=\s*(.+)", block)
        rows.append({
            "GSM":   gsm_m.group(1).strip() if gsm_m else "",
            "title": title_m.group(1).strip() if title_m else "",
            "characteristics": " | ".join(c.strip() for c in char_m[:3])
        })
    result["samples"] = rows
    return result


def build_context(gse, sra_rows, soft_info, qc_summary, max_chars=30000):
    """组装 LLM 上下文字符串"""
    parts = [f"数据集: {gse}"]

    if soft_info:
        parts.append(f"\n标题: {soft_info.get('title', 'N/A')}")
        parts.append(f"\n研究摘要:\n{soft_info.get('summary', 'N/A')[:3000]}")
        parts.append(f"\n整体设计:\n{soft_info.get('overall_design', 'N/A')[:2000]}")

    if sra_rows:
        total = len(sra_rows)
        strats = {}
        layouts = {}
        for r in sra_rows:
            s = r.get("LibraryStrategy", "?")
            strats[s] = strats.get(s, 0) + 1
            l_ = r.get("LibraryLayout", "?")
            layouts[l_] = layouts.get(l_, 0) + 1
        total_spots = sum(int(r.get("spots") or 0) for r in sra_rows)
        total_gb = sum(int(r.get("bases") or 0) for r in sra_rows) / 1e9
        parts.append(
            f"\nSRA 信息: {total} 个 Run | "
            f"策略: {', '.join(f'{k}:{v}' for k,v in strats.items())} | "
            f"文库布局: {', '.join(f'{k}:{v}' for k,v in layouts.items())} | "
            f"总 spots: {total_spots:,} | 总 bases: {total_gb:.1f} GB"
        )
        # 样本表前 50 行
        sample_lines = ["Run\tGSM\tTitle\tStrategy\tLayout\tspots"]
        for r in sra_rows[:50]:
            sample_lines.append(
                f"{r.get('Run','')}\t{r.get('GSM','')}\t{r.get('Title','')[:60]}\t"
                f"{r.get('LibraryStrategy','')}\t{r.get('LibraryLayout','')}\t{r.get('spots','')}"
            )
        parts.append("\n--- 样本表（前50行）---\n" + "\n".join(sample_lines))

    if soft_info and soft_info.get("samples"):
        char_lines = ["GSM\t标题\t样本特征"]
        for s in soft_info["samples"][:50]:
            char_lines.append(f"{s['GSM']}\t{s['title'][:50]}\t{s['characteristics'][:80]}")
        parts.append("\n--- GEO 样本特征（前50行）---\n" + "\n".join(char_lines))

    if qc_summary:
        parts.append(f"\n--- QC 快照 ---\n{qc_summary}")

    context = "\n".join(parts)
    if len(context) > max_chars:
        context = context[:max_chars] + "\n...[上下文已截断]"
    return context


SYSTEM_PROMPT = """你是资深生信数据审阅专家。请基于提供的 GEO/SRA 元数据，
输出中文数据解读卡，必须严格包含以下 JSON 字段（返回合法 JSON，不加 markdown 代码块）：
{
  "gse": "GSE编号",
  "purpose": "研究目的（1-3句）",
  "design": "实验设计简述（分组、时间点、干预等）",
  "groups": ["组1描述", "组2描述"],
  "n_samples": 样本数,
  "platform": "测序平台",
  "sample_source_type": "样本来源类型（如: 原代肿瘤组织/细胞系/类器官/PDX/血液/PBMC/正常组织）",
  "tissue_or_cell_line": "具体组织部位或细胞系名称（如: 人GBM组织、T98G细胞系、星形胶质细胞 等）",
  "tumor_match_score": 0-100整数（与研究者指定肿瘤类型的匹配度；若无背景提示则评估数据通用性）,
  "tumor_match_reason": "匹配度评分理由（说明样本是原发组织、细胞系、PDX还是其他，是否与指定肿瘤相关）",
  "quality_risks": ["风险1", "风险2"],
  "batch_effect_notes": "批次效应提示",
  "reusability_score": 0-100整数（综合可复用性：样本质量+数据完整性+分组合理性）,
  "reusability_reason": "综合评分理由"
}"""


def call_llm(context, gse, llm_cfg, api_key):
    """调用 LLM API，返回 JSON 字符串或 None"""
    if not HAS_OPENAI:
        return None, "openai 未安装，请 pip install openai"

    client = OpenAI(
        api_key=api_key,
        base_url=llm_cfg.get("api_base", "https://api.openai.com/v1"),
    )
    proxy = llm_cfg.get("proxy", "")
    if proxy:
        import httpx
        client._client = httpx.Client(proxies={"https://": proxy})

    user_msg = f"请分析以下数据集：\n\n{context}"
    try:
        response = client.chat.completions.create(
            model=llm_cfg.get("model", "gpt-4o-mini"),
            messages=[
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user",   "content": user_msg},
            ],
            temperature=float(llm_cfg.get("temperature", 0.2)),
            timeout=float(llm_cfg.get("timeout_sec", 120)),
        )
        return response.choices[0].message.content, None
    except Exception as e:
        return None, str(e)


def parse_llm_response(raw_text, gse):
    """从 LLM 返回文本中提取 JSON，兜底返回最小结构"""
    if not raw_text:
        return {"gse": gse, "error": "empty response"}
    # 去掉可能的 ```json ... ``` 包裹
    text = re.sub(r"```json\s*", "", raw_text)
    text = re.sub(r"```\s*", "", text)
    # 尝试解析第一个 { ... }
    m = re.search(r"\{.*\}", text, re.DOTALL)
    if m:
        try:
            return json.loads(m.group(0))
        except Exception:
            pass
    return {"gse": gse, "raw_response": raw_text[:2000], "parse_error": True}


def write_outputs(parsed, gse, output_dir, context_used=""):
    """写出 .json 和 .md 解读卡"""
    os.makedirs(output_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # JSON
    parsed["_generated_at"] = ts
    json_path = os.path.join(output_dir, f"{gse}.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(parsed, f, ensure_ascii=False, indent=2)

    # Markdown 解读卡
    md_path = os.path.join(output_dir, f"{gse}.md")
    score = parsed.get("reusability_score", "N/A")
    score_bar = "🟢" if isinstance(score, int) and score >= 70 else \
                "🟡" if isinstance(score, int) and score >= 40 else "🔴"
    tumor_match = parsed.get("tumor_match_score", "N/A")
    tumor_bar = "🟢" if isinstance(tumor_match, int) and tumor_match >= 70 else \
                "🟡" if isinstance(tumor_match, int) and tumor_match >= 40 else "🔴"

    md = f"""# {gse} 数据集解读卡

> 生成时间: {ts} | 可复用性: {score_bar} **{score}/100** | 瘤种匹配度: {tumor_bar} **{tumor_match}/100**

## 研究目的
{parsed.get('purpose', 'N/A')}

## 实验设计
{parsed.get('design', 'N/A')}

## 分组结构
{chr(10).join(f'- {g}' for g in parsed.get('groups', []))}

## 基本信息
- 样本数: {parsed.get('n_samples', 'N/A')}
- 测序平台: {parsed.get('platform', 'N/A')}

## 样本来源与瘤种匹配
- **样本类型**: {parsed.get('sample_source_type', 'N/A')}
- **具体组织/细胞系**: {parsed.get('tissue_or_cell_line', 'N/A')}
- **匹配度评分**: {tumor_bar} **{tumor_match}/100**
- **评分理由**: {parsed.get('tumor_match_reason', 'N/A')}

## 质量与批次效应风险
{chr(10).join(f'⚠️ {r}' for r in parsed.get('quality_risks', []))}

**批次效应提示**: {parsed.get('batch_effect_notes', 'N/A')}

## 可复用性评分
**{score}/100** — {parsed.get('reusability_reason', 'N/A')}
"""
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(md)

    print(f"[LLM] 解读卡已写出: {json_path}")
    print(f"[LLM] Markdown 已写出: {md_path}")
    return json_path, md_path


# ──────────────────────────────────────────────
# 主流程
# ──────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="GEO 数据集 LLM 智能解读")
    p.add_argument("--gse",        required=True, help="GSE 编号，如 GSE242225")
    p.add_argument("--project",    default="",    help="项目名（用于 output_dir 路径拼接）")
    p.add_argument("--sra_info",   default="workflow/resources/homo/",
                   help="SraRunInfo_*.csv 所在目录")
    p.add_argument("--output_dir", default=None,  help="解读卡输出目录（默认: result/{project}/00_data_intel）")
    p.add_argument("--qc_dir",     default=None,  help="03_Align_Filter 目录（提取 QC 摘要）")
    p.add_argument("--species",    default="homo", help="物种（用于 qc_dir 路径）")
    p.add_argument("--force",      action="store_true", help="强制重新生成（覆盖缓存）")
    p.add_argument("--llm_config", default="config/llm.yaml", help="llm.yaml 路径")
    p.add_argument("--env_file",   default="config/.env",     help=".env 密钥文件路径")
    return p.parse_args()


def main():
    args = parse_args()
    gse = args.gse.strip().upper()

    # 输出目录
    if args.output_dir:
        output_dir = args.output_dir
    elif args.project:
        output_dir = os.path.join("result", args.project, "00_data_intel")
    else:
        output_dir = "result/00_data_intel"

    # 缓存检查
    json_path = os.path.join(output_dir, f"{gse}.json")
    if not args.force and os.path.exists(json_path):
        print(f"[LLM] {gse} 解读卡已存在（缓存），跳过（--force 可覆盖）")
        return

    # 加载配置
    load_env(args.env_file)
    llm_cfg = load_llm_config(args.llm_config)
    api_key = os.environ.get("LLM_API_KEY", "")
    if not api_key:
        print("[ERROR] 未找到 LLM_API_KEY，请在 config/.env 中设置")
        sys.exit(1)

    # 加载数据
    print(f"[LLM] 处理 {gse}...")
    sra_rows = load_sra_info(args.sra_info, gse)
    print(f"[LLM] 读取 SRA 信息: {len(sra_rows)} 条")

    print(f"[LLM] 下载 GEO SOFT 元数据...")
    soft_text = fetch_geo_soft_summary(gse)
    soft_info = parse_soft(soft_text,
                           max_sample_rows=llm_cfg.get("max_context_chars", 50))

    qc_summary = None
    if args.qc_dir and os.path.isdir(args.qc_dir):
        qc_summary = load_qc_summary(args.qc_dir, gse, args.species)

    # 组装上下文
    max_chars = int(llm_cfg.get("max_context_chars", 30000))
    context = build_context(gse, sra_rows, soft_info, qc_summary, max_chars)
    print(f"[LLM] 上下文长度: {len(context)} 字符")

    # 调用 API
    print(f"[LLM] 调用 {llm_cfg.get('model', '未配置')}...")
    t0 = time.time()
    raw, err = call_llm(context, gse, llm_cfg, api_key)
    elapsed = time.time() - t0

    if err:
        print(f"[LLM] API 调用失败: {err}")
        # 降级：写出失败记录
        fallback = {"gse": gse, "error": err, "offline_mode": True,
                    "_generated_at": datetime.now().isoformat()}
        os.makedirs(output_dir, exist_ok=True)
        with open(json_path, "w") as f:
            json.dump(fallback, f, ensure_ascii=False, indent=2)
        print(f"[LLM] 失败记录已写出: {json_path}")
        return

    print(f"[LLM] API 响应 ({elapsed:.1f}s)，解析结果...")
    parsed = parse_llm_response(raw, gse)
    write_outputs(parsed, gse, output_dir)


if __name__ == "__main__":
    main()
