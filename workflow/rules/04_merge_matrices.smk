"""
rule 04_merge_matrices
=======================
功能: 扫描所有 GSE 的 StringTie 定量结果（gene_abund.tab），
      按物种和注释类型分组，合并为表达量矩阵 CSV（TPM），
      支持可选的 QC 过滤（仅保留 Passed 样本）。

依赖:
  - rule filter_alignment 完成后才能执行
  - 依赖脚本: workflow/scripts/04_merge_matrices.py

输   入: result/{project}/03_Align_Filter/alignment_quality.csv（含 filter 信息）
输   出: result/{project}/04_merge_matrices/Merge_finished.txt
         result/{project}/04_merge_matrices/Matrices_<timestamp>/
             human_mRNA_genecode_stringtie_matrix.csv
             human_lncRNA_GENCODE_stringtie_matrix.csv
             human_lncRNA_NONCODE_stringtie_matrix.csv
             human_eRNA_EnhancerAtlas_stringtie_matrix.csv
             human_eRNA_Ensembl_stringtie_matrix.csv
             human_eRNA_FANTOM5_stringtie_matrix.csv
             human_miRNA_miRBase_stringtie_matrix.csv
             human_miRNA_MirGeneDB_stringtie_matrix.csv

脚本说明 (04_merge_matrices.py):
  - 递归搜索 gene_abund.tab 文件
  - 根据路径自动识别物种（homo→human, mouse→mouse）
  - 根据路径匹配注释类型（mRNA/eRNA/lncRNA/miRNA × 不同数据库）
  - 输出目录加时间戳避免覆盖历史结果
"""

import os

def _get_align_filter_base(wildcards):
    """
    返回 03_Align_Filter 目录（不含物种子目录），
    让 04_merge_matrices.py 能从相对路径中检测物种名（homo/mouse）。
    """
    return os.path.join("result", wildcards.project, "03_Align_Filter")

rule merge_matrices:
    input:
        script       = "workflow/scripts/04_merge_matrices.py",
        filter_csv   = "result/{project}/03_Align_Filter/alignment_quality.csv",
        filter_marker = "result/{project}/03_Align_Filter/Filter_finished.txt"
    output:
        marker = "result/{project}/04_merge_matrices/Merge_finished.txt"
    params:
        inputdir  = _get_align_filter_base,
        outputdir = "result/{project}/04_merge_matrices/Matrices",
        gtf_base  = config.get("gtf_base", "/home/public_software_annotation")
    log:
        "logs/{project}/04_merge_matrices.log"
    threads: config.get("merge_threads", 1)
    shell:
        """
        mkdir -p result/{wildcards.project}/04_merge_matrices
        python {input.script} \
            --inputdir  {params.inputdir} \
            --outputdir {params.outputdir} \
            --filter_csv {input.filter_csv} \
            --gtf_base  {params.gtf_base} > {log} 2>&1

        # 04 脚本将 Merge_finished.txt 写到 outputdir 的父目录（即 04_merge_matrices/）
        # 若文件未在预期位置，手动生成
        if [ ! -f {output.marker} ]; then
            echo "Merge_finished.txt not found at expected path, creating." >> {log}
            touch {output.marker}
        fi
        """
