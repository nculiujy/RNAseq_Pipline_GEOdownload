# RNAseq_GEO workflow/scripts — 按执行顺序整理

GSE 原始数据 → 各批次（GSE）样本表达量矩阵 的脚本序列。
对应 Yuseq_Pipline 的 `result/3_Align_Filter/` 批次产出（每个 GSE 一个批次目录 + 每日期一个合并矩阵目录）。

## 执行顺序

| 序号 | 脚本 | 输入 | 输出 |
|------|------|------|------|
| 00 | `00_gse_to_srr.py` | GSE 列表 CSV（任意含 `GSE_ID` 列） | `SRR_Acc_List_<GSE>.txt` + `SraRunInfo_<GSE>.csv` + `GSE_SRR_summary.csv` + `ALL_rnaseq_SRR.txt` |
| 01 | `01_download_sra.py` | 含 `SRR.txt` 的目录树（递归扫描） | `.sra` 文件 → `result_dir/<species>/rawdata/<GSE>/`（**只下载不解压**，.sra 由 02 解压） |
| 02 | `02_dataset_pipeline.py` | `.sra`（每 GSE 一个任务） | fastq→clean→bam→stringtie 定量（8 类注释）+ `dataset_finished.txt` |
| 03 | `03_filter_alignment.py` | `3_Align_Filter/` 下 `QC_results.log` | `alignment_quality.csv`（比对率 ≥70% 判 Passed） |
| 04 | `04_merge_matrices.py` | 定量结果目录（+ 可选 filter csv） | 每类注释一张表达量矩阵（`*_stringtie_matrix.csv`） |

## 数据流

```
GSE 列表 CSV
   │ 00_gse_to_srr.py
   ▼
SRR_Acc_List_<GSE>.txt ──复制/链接为──▶ <rawdata>/<species>/rawdata/<GSE>/SRR.txt
   │ 01_download_sra.py (prefetch，只下载)
   ▼
result/1_download_geo/<species>/rawdata/<GSE>/*.sra
   │ 02_dataset_pipeline.py (fasterq-dump→fastp→HISAT2→StringTie，每 GSE 一批)
   ▼
result/3_Align_Filter/<species>/<GSE>/hisat2file/<SRR>/  +  dataset_finished.txt
   │ 03_filter_alignment.py（可选，产出比对率过滤表）
   ▼
result/3_Align_Filter/<species>/<GSE>/alignment_quality.csv
   │ 04_merge_matrices.py（合并全部 GSE → 每日期目录）
   ▼
result/3_Align_Filter/<YYYY.MM.DD>/<rna>_<anno>_stringtie_matrix.csv   ← 批次表达量矩阵
```

## Snakemake 映射建议（rules/ 下新建）

```python
# rule 00: 生成 SRR 列表（可手动执行一次，或用 --csv 批量）
# rule 01: download_sra    —— 通配符 {species}/{gse}，输出 1_download_geo/finished.txt
# rule 02: per_dataset     —— 通配符 {species}/{gse}，输出 3_Align_Filter/{species}/{gse}/dataset_finished.txt
# rule 03: filter_align    —— 依赖 rule02 输出，产出 alignment_quality.csv
# rule 04: merge           —— expand 所有 {species}/{gse} 的 dataset_finished.txt，产出 Merge_finished.txt + Matrices/
```

参考现有实现：`Yuseq_Pipline/workflow/rules/3_Align_Filter.smk`（per_dataset_pipeline + Merge 两个 rule）。

## 关键参数

- 02_dataset_pipeline.py：`--threads`（每样本线程，默认8）、`--maxparallel`（数据集内并行样本数，默认4）
- HISAT2 索引 / GTF 注释基目录：脚本内 `GTF_BASE=/home/public_software_annotation`，species=homo/mouse
- 注释类型（8 类）：mRNA(GENCODE)、lncRNA(GENCODE/NONCODE)、eRNA(EnhancerAtlas/Ensembl/FANTOM5)、miRNA(miRBase/MirGeneDB)

## 已知流程要点（2026-08-08 更新）

- ✅ 流程覆盖 RNA-seq 基础环节：下载→解压→QC/trim→HISAT2 比对（--dta + RG）→去重→StringTie 定量（-e -B，8类注释）→比对率过滤→表达矩阵
- ✅ 双端/单端自动识别；断点续跑（fastq/bam/gtf 存在即跳过）
- ✅ **P0 已修复**：01 原自带 fasterq-dump 解压并删除 .sra，与 02 的输入冲突；现 01 只下载、02 负责解压
- ⚠️ P1-1 去重策略：MarkDuplicates REMOVE_DUPLICATES=true 对 RNA-seq 有争议（可能损失真实信号），建议按研究设计确认
- ✅ **P1-2 已修复**：02 有样本失败时现以 `sys.exit(1)` 非零退出，Snakemake 正确标记该 GSE 任务失败；断点续跑用 `--rerun-incomplete`
- ✅ **P2-1 已修复**：02 新增 `--strandedness` 参数（`unstranded`/`forward`/`reverse`），透传 HISAT2 `--rna-strandness` 和 StringTie `--rf`/`--fr`；通过 `config.yaml` 全局配置
- ⚠️ P2-2 fastp 建议：默认参数可跑，正式分析建议显式配置 adapter/质量阈值；如需 FastQC 报告可另加
- ✅ **P3 已修复**：fasterq-dump 已改用 `--split-3`（比 `--split-files` 更稳，兼容异常 paired 记录）

## archive/（已弃用，保留备查）

| 旧文件 | 被取代原因 |
|--------|-----------|
| `2_QC_fqfile.pl` | fastp QC 已并入 02_dataset_pipeline.py（snakefile 注释：2_QC 已合并） |
| `3_1_Align.py` | 全局批处理版 HISAT2，被 02 内置 run_hisat2 取代 |
| `3_2_Quant.py` | 全局批处理版 StringTie，被 02 内置 run_stringtie 取代 |
| `3_4_merge.py` | 旧版合并（filter 必填），被 04_merge_matrices.py（filter 可选+自动物种检测）取代 |
