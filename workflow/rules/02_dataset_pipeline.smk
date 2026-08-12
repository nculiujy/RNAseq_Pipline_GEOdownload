"""
rule 02_dataset_pipeline（per-GSE，批次串行化 + 即跑即删）
============================================================
功能: 对单个 GSE 数据集执行完整的分析流水线：
        fasterq-dump 解压 → fastp QC → HISAT2 比对 →
        Picard MarkDuplicates → StringTie 定量（8 类注释）→
        ★ 删除 .sra 和 fastq（即跑即删，释放磁盘）→
        写出 dataset_finished.txt

核心改动（对应批次续跑方案 v2）：
  - input.download_marker = per-GSE 的 .download_done（而非全局 download_finished.txt）
  - resources: gse_slots=1（同一时刻只处理一个 GSE，与 download_gse 串行化）
  - shell 末尾删除 .sra（+ 可选 bam）——只有 02 成功才执行到这里

通配符:
  {project}  = config 中的 project_name（如 GBM_homo）
  {species}  = 物种（homo | mouse）
  {gse}      = GSE 编号（如 GSE119834）

运行方式:
  snakemake -j 32 --resources gse_slots=1
"""

import os

def _get_proj_02(wildcards):
    for p in config["projects"]:
        if p["project_name"] == wildcards.project.split("/")[0]:
            return p
    return {}

def get_sra_dir(wildcards):
    return os.path.join(
        "result", wildcards.project,
        "01_download_sra", wildcards.species, "rawdata", wildcards.gse
    )

# ──────────────────────────────────────────────
# 核心 rule：每个 {project}/{species}/{gse} 一个任务
# ──────────────────────────────────────────────
rule dataset_pipeline:
    input:
        script          = "workflow/scripts/02_dataset_pipeline.py",
        # ★ 依赖 per-GSE 的 .download_done（不再依赖全局 download_finished.txt）
        download_marker = "result/{project}/01_download_sra/{species}/rawdata/{gse}/.download_done"
    output:
        marker = "result/{project}/03_Align_Filter/{species}/{gse}/dataset_finished.txt"
    params:
        dataset_id  = lambda wildcards: wildcards.gse,
        species     = lambda wildcards: wildcards.species,
        sra_dir     = lambda wildcards: get_sra_dir(wildcards),
        fastq_dir   = lambda wildcards: os.path.join(
                          "result", wildcards.project,
                          "01_download_sra", wildcards.species, "rawdata", wildcards.gse),
        qc_dir      = lambda wildcards: os.path.join(
                          "result", wildcards.project,
                          "02_QC", wildcards.species, wildcards.gse),
        align_dir   = lambda wildcards: os.path.join(
                          "result", wildcards.project,
                          "03_Align_Filter", wildcards.species, wildcards.gse),
        threads           = config.get("pipeline_threads",   8),
        maxparallel       = config.get("pipeline_parallel",  4),
        strandedness      = config.get("strandedness", "unstranded"),
        picard_jar        = config.get("picard_jar", "workflow/env/picard-2.18.2/picard.jar"),
        keep_bam          = config.get("cleanup", {}).get("keep_bam", False),
        keep_sra          = config.get("cleanup", {}).get("keep_sra", False),
        sample_chunk_size = config.get("sample_chunk_size", 0),
        anno_base         = config.get("anno_base", ""),
        hisat2_index      = lambda wildcards: config.get("hisat2_index", {}).get(wildcards.species, ""),
        remove_duplicates = "true" if config.get("remove_duplicates", True) else "false",
    log:
        "logs/{project}_{species}_{gse}_02_dataset_pipeline.log"
    threads: config.get("pipeline_threads", 8) * config.get("pipeline_parallel", 4)
    resources:
        gse_slots = 1    # ★ 串行化：同一时刻只处理一个 GSE
    shell:
        """
        mkdir -p {params.qc_dir} {params.align_dir}

        # 将 cleanup 策略通过环境变量传入脚本（滚动窗口模式下 cleanup_chunk 使用）
        export KEEP_SRA={params.keep_sra}
        export KEEP_BAM={params.keep_bam}

        python {input.script} \
            --dataset_id         {params.dataset_id} \
            --species            {params.species} \
            --sra_dir            {params.sra_dir} \
            --fastq_dir          {params.fastq_dir} \
            --qc_dir             {params.qc_dir} \
            --align_dir          {params.align_dir} \
            --threads            {params.threads} \
            --maxparallel        {params.maxparallel} \
            --strandedness       {params.strandedness} \
            --picard_jar         "{params.picard_jar}" \
            --sample_chunk_size  {params.sample_chunk_size} \
            --anno_base          "{params.anno_base}" \
            --hisat2_index       "{params.hisat2_index}" \
            --remove_duplicates  {params.remove_duplicates} \
            --output_marker      {output.marker} > {log} 2>&1

        # ★ 即跑即删（shell 层兜底清理）：
        #   - 分块模式(sample_chunk_size>0)：脚本内 cleanup_chunk 已逐批清理，此处兜底
        #   - 非分块模式(sample_chunk_size=0)：全部由此处清理
        if [ "{params.keep_sra}" != "True" ]; then
            echo "[Cleanup] 删除残余 .sra / fastq 文件..." >> {log}
            rm -f {params.sra_dir}/*.sra
            rm -f {params.fastq_dir}/*.fastq {params.fastq_dir}/*.fastq.gz
            rm -f {params.qc_dir}/*.clean.fastq
        fi

        # 可选：删除 bam（保留表达矩阵输入，节省磁盘）
        if [ "{params.keep_bam}" != "True" ]; then
            echo "[Cleanup] 删除残余 bam 文件..." >> {log}
            rm -rf {params.align_dir}/hisat2file
        fi
        """
