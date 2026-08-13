"""
rule download_gse（per-GSE，批次串行化）
=========================================
功能: 对单个 GSE 运行下载 .sra 文件。
      优先使用 aria2c 多连接加速（比 prefetch 快 5-10 倍），
      若 aria2c 不可用则自动回退到 prefetch。
      通过 resources: gse_slots=1 实现批次串行化——
      同一时刻只有一个 GSE 在下载，避免磁盘/带宽峰值过高。

依赖脚本: workflow/scripts/01_download_sra.py（支持 --gse 参数）
通配符:
  {project} — config 中的 project_name
  {species}  — 物种
  {gse}      — GSE 编号（如 GSE119834）
输   出: result/{project}/01_download_sra/{species}/rawdata/{gse}/.download_done

运行方式:
  snakemake -j 32 --resources gse_slots=1
  （串行：同一时刻只有 1 个 GSE 在下载；调大 gse_slots 可改为并行几个）

aria2c 安装（若未安装）:
  conda install -c conda-forge aria2
"""

import os

def _get_proj_per_gse(wildcards):
    for p in config["projects"]:
        if p["project_name"] == wildcards.project.split("/")[0]:
            return p
    return {}

def get_rawdata_dir_pergse(wildcards):
    p = _get_proj_per_gse(wildcards)
    return p.get("rawdata_dir", f"workflow/resources/{wildcards.species}")

rule download_gse:
    input:
        script      = "workflow/scripts/01_download_sra.py",
        rawdata_dir = get_rawdata_dir_pergse   # 目录作为输入，存在即可
    output:
        marker = "result/{project}/01_download_sra/{species}/rawdata/{gse}/.download_done"
    params:
        rawdata_dir  = get_rawdata_dir_pergse,
        result_dir   = lambda wildcards: os.path.join(
            "result", wildcards.project,
            "01_download_sra", wildcards.species, "rawdata"
        ),
        gse          = lambda wildcards: wildcards.gse,
        # aria2c 参数从 config 读取，提供默认值
        connections  = config.get("aria2c_connections", 16),
        split        = config.get("aria2c_split",       16),
        chunk        = config.get("aria2c_chunk",       "1M"),
        timeout      = config.get("aria2c_timeout",     300),
    log:
        "logs/{project}/download_{species}_{gse}.log"
    threads: config.get("download_threads", 8)
    resources:
        gse_slots = 1    # ★ 串行化：同一时刻只下载一个 GSE
    shell:
        """
        mkdir -p {params.result_dir}
        export DOWNLOAD_THREADS={threads}
        export ARIA2_CONNECTIONS={params.connections}
        export ARIA2_SPLIT={params.split}
        export ARIA2_CHUNK={params.chunk}
        export ARIA2_TIMEOUT={params.timeout}
        python {input.script} \
            {params.rawdata_dir} \
            {params.result_dir} \
            --gse {params.gse} > {log} 2>&1

        # 无论成功还是失败都写 .download_done（失败的 SRR 由 02 的跳过逻辑处理）
        touch {output.marker}
        echo "[Download] {params.gse} 完成（含部分失败时仍继续）" >> {log}
        """
