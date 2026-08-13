"""
rule 03_filter_alignment
=========================
功能: 扫描所有 GSE 子目录下的 QC_results.log（HISAT2 比对日志），
      提取比对率，按阈值（默认 70%）标记 Passed/Failed，
      输出汇总 CSV 文件。

依赖:
  - rule dataset_pipeline 完成后才能执行
  - 依赖脚本: workflow/scripts/03_filter_alignment.py

输   入: 所有 GSE 的 dataset_finished.txt（expand）
输   出: result/{project}/03_Align_Filter/alignment_quality.csv
         result/{project}/03_Align_Filter/Filter_finished.txt

脚本说明 (03_filter_alignment.py):
  - 递归扫描 --inputdir 下所有 QC_results.log
  - 提取最后一个百分比数值作为比对率
  - 写出 alignment_quality.csv（Sample_ID, Alignment_Rate, Passed, Path）
"""

import os
import csv as _csv

def _pname_from_wildcard(project_wildcard):
    """
    从 {project} wildcard 提取真正的 project_name。
    当 run_id 被合并入 project 时，格式为 {pname}/{run_id}，
    project_name 是第一个路径组件。
    """
    return project_wildcard.split("/")[0]


def _get_gse_dataset_markers(wildcards):
    """
    收集该 project 下所有 GSE 的 dataset_finished.txt 路径，
    作为 filter_alignment rule 的聚合输入，确保所有数据集都处理完成。
    """
    # 从 Snakefile 的 TARGET_FILES 中反向查找本 project 对应的 dataset_finished.txt
    # 直接扫描已有的 dataset_finished.txt 文件更可靠
    import glob as _glob
    base = os.path.join("result", wildcards.project, "03_Align_Filter")
    markers = _glob.glob(os.path.join(base, "*", "*", "dataset_finished.txt"))
    # 若没有已完成文件，则从 Snakefile TARGET_FILES 中找
    if not markers:
        pname   = _pname_from_wildcard(wildcards.project)
        for proj in config["projects"]:
            if proj["project_name"] != pname:
                continue
            species   = proj["species"]
            rawdata   = proj.get("rawdata_dir", f"workflow/resources/{species}")
            # 使用与 Snakefile 相同的 GSE 列表读取逻辑
            batch_f   = config.get("batch_file", None)
            if batch_f and os.path.exists(batch_f):
                table = batch_f
            elif os.path.exists(os.path.join(rawdata, "batch_input.txt")):
                table = os.path.join(rawdata, "batch_input.txt")
            else:
                table = os.path.join(rawdata, "SRR_table.txt")
            gse_list = []
            if os.path.exists(table):
                import re as _re_smk
                with open(table) as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        m = _re_smk.search(r"(GSE\d+)", line, _re_smk.I)
                        if m:
                            gse = m.group(1).upper()
                            if gse not in gse_list:
                                gse_list.append(gse)
            for gse in gse_list:
                markers.append(
                    os.path.join(
                        "result", wildcards.project,
                        "03_Align_Filter", species, gse, "dataset_finished.txt"
                    )
                )
    return markers


def _get_align_filter_dir(wildcards):
    pname = _pname_from_wildcard(wildcards.project)
    for proj in config["projects"]:
        if proj["project_name"] == pname:
            return os.path.join(
                "result", wildcards.project,
                "03_Align_Filter", proj["species"]
            )
    return os.path.join("result", wildcards.project, "03_Align_Filter")

rule filter_alignment:
    input:
        script           = "workflow/scripts/03_filter_alignment.py",
        dataset_markers  = _get_gse_dataset_markers   # 聚合：等待所有 GSE 完成
    output:
        csv    = "result/{project}/03_Align_Filter/alignment_quality.csv",
        marker = "result/{project}/03_Align_Filter/Filter_finished.txt"
    params:
        inputdir  = _get_align_filter_dir,
        outputdir = "result/{project}/03_Align_Filter",
        cutoff    = config.get("alignment_rate_cutoff", 70.0)
    log:
        "logs/{project}/03_filter_alignment.log"
    threads: config.get("filter_threads", 1)
    shell:
        """
        mkdir -p {params.outputdir}
        python {input.script} \
            --inputdir  {params.inputdir} \
            --outputdir {params.outputdir} \
            --cutoff    {params.cutoff} > {log} 2>&1

        # 03 脚本自己会写 Filter_finished.txt；若不在预期位置则移动
        if [ ! -f {output.marker} ]; then
            if [ -f {params.outputdir}/Filter_finished.txt ]; then
                mv {params.outputdir}/Filter_finished.txt {output.marker}
            else
                touch {output.marker}
            fi
        fi
        """
