# RNAseq_GEO — GEO 公共数据集 RNA-seq 全自动分析流程

基于 **Snakemake** 的批量 RNA-seq 流程：GSE 号列表 → SRR 获取 → 下载 → QC → HISAT2 比对 → StringTie 多注释定量 → 表达矩阵，配套 **Streamlit UI** 完成全部配置与监控。

---

## 1. conda 环境依赖

```bash
conda env create -f environment.yml
conda activate RNAseq_Pipline
```

环境包含：Python 3.13、hisat2 2.2.2、fastp 1.1.0、stringtie 3.0.3、samtools 1.21、sra-tools 3.2.1、snakemake 9.17.2、streamlit 1.59.1、pandas / plotly / openai / duckdb 等。

> 国内网络建议先配置 conda 镜像源（清华 TUNA）。
> sra-tools 需 ≥3.x（2.x 有 TLS 兼容问题无法连接 NCBI）。

## 2. 注释文件放置位置

以下文件 conda 装不了，需自行放到指定路径：

```bash
# Picard jar
mkdir -p workflow/env/picard-2.18.2
cp <你的路径>/picard.jar workflow/env/picard-2.18.2/

# 参考基因组索引 + GTF 注释（hg38）
workflow/anno/homo/
├── humanHisat2Index/GRCh38/     # HISAT2 索引（.ht2 文件）
├── mRNAanno/                    # gencode.v44.mRNA.annotation.gtf
└── ncRNAanno/                   # eRNA/lncRNA/miRNA 共 7 个 GTF
```

## 3. UI 启动

```bash
conda activate RNAseq_Pipline
streamlit run streamlit_app/app.py --server.port 8501
```

浏览器打开 `http://localhost:8501` 后，项目配置（GSE 列表、批次、并发、链特异性、LLM 等）均可在 UI 内完成，无需手动改配置文件。

---

## 项目结构

```
RNAseq_GEO/
├── Snakefile                # Snakemake 主入口（00-04 模块）
├── config/                  # 全局配置（config.yaml / llm.yaml / .env）
├── workflow/
│   ├── scripts/             # 步骤脚本（00-10）
│   ├── rules/               # Snakemake rule
│   ├── resources/{species}/ # GSE 列表 / 批次文件
│   ├── anno/                # 参考注释（自行放置）
│   └── env/                 # Picard jar（自行放置）
├── streamlit_app/           # UI 应用
├── run_all.sh               # 批次自动续跑启动脚本
├── environment.yml          # conda 环境
└── result/{project}/{run_id}/  # 运行产物（不进版本库）
```
