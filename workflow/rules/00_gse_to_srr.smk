"""
rule 00_gse_to_srr
==================
功能: 读取含 GSE_ID 列的 CSV，批量拉取每个 GSE 的 SRR Acc List，
      产出 SRR_Acc_List_<GSE>_rnaseq.txt + SraRunInfo_<GSE>.csv +
      更新后的 GSE_SRR_summary.csv + ALL_rnaseq_SRR.txt。

依赖脚本: workflow/scripts/00_gse_to_srr.py
输   出: result/{project}/00_gse_to_srr/GSE_SRR_summary_updated.csv

注意:
  - 此步骤一般手动执行一次，后续重跑会自动跳过已有结果（幂等）
  - 若不想通过 Snakemake 执行，可直接运行：
      python workflow/scripts/00_gse_to_srr.py \\
          --csv workflow/resources/homo/GSE_SRR_summary.csv \\
          -o workflow/resources/homo/
"""

import os

def _get_proj(wildcards):
    for p in config["projects"]:
        if p["project_name"] == wildcards.project:
            return p
    return {}

def get_gse_list_csv(wildcards):
    p = _get_proj(wildcards)
    return p.get("gse_list_csv", "")

rule gse_to_srr:
    input:
        script = "workflow/scripts/00_gse_to_srr.py",
        gse_csv = get_gse_list_csv
    output:
        summary = "result/{project}/00_gse_to_srr/GSE_SRR_summary_updated.csv"
    params:
        outdir = "result/{project}/00_gse_to_srr"
    log:
        "logs/{project}/00_gse_to_srr.log"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        python {input.script} \
            --csv {input.gse_csv} \
            -o {params.outdir} > {log} 2>&1

        # 将更新后的汇总文件重命名为输出目标（方便 Snakemake 追踪）
        if [ -f {params.outdir}/GSE_SRR_summary.csv ]; then
            cp {params.outdir}/GSE_SRR_summary.csv {output.summary}
        else
            echo "WARNING: GSE_SRR_summary.csv not found in {params.outdir}" >> {log}
            touch {output.summary}
        fi
        """
