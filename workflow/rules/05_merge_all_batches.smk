"""
rule 05_merge_all_batches — 跨批次矩阵整合
=============================================
功能: 扫描 result/{project}/*/03_Align_Filter/{species}/ 下所有 run_id 的
      gene_abund.tab，合并为全量表达矩阵 + sample_manifest.csv。

输出:
  result/{project}/00_final_matrices/Merge_all_finished.txt
  result/{project}/00_final_matrices/Matrices_all/*_matrix.csv
  result/{project}/00_final_matrices/sample_manifest.csv

使用方式:
  - 所有批次跑完后手动触发: snakemake result/GBM_homo/00_final_matrices/Merge_all_finished.txt
  - 或在 Snakefile TARGET_FILES 中追加此目标
  - 重复执行幂等（覆盖旧文件）

注意:
  此 rule 不依赖具体 run_id 的 marker（因为需要跨 run_id 扫描），
  建议在确认所有批次完成后手动触发。
"""

import os


def _get_species_for_project(wildcards):
    """从 config 中获取 project 对应的 species"""
    pname = wildcards.project.split("/")[0]
    for proj in config["projects"]:
        if proj["project_name"] == pname:
            return proj["species"]
    return "homo"  # fallback


rule merge_all_batches:
    input:
        script = "workflow/scripts/05_merge_all_batches.py"
    output:
        marker = "result/{project}/00_final_matrices/Merge_all_finished.txt"
    params:
        result_base = lambda wildcards: os.path.join("result", wildcards.project),
        species     = _get_species_for_project,
        output_dir  = lambda wildcards: os.path.join("result", wildcards.project, "00_final_matrices"),
        cutoff      = config.get("alignment_rate_cutoff", 70.0)
    log:
        "logs/{project}/05_merge_all_batches.log"
    threads: config.get("merge_threads", 1)
    shell:
        """
        mkdir -p {params.output_dir}
        python {input.script} \
            --result_base {params.result_base} \
            --species     {params.species} \
            --output_dir  {params.output_dir} \
            --cutoff      {params.cutoff} > {log} 2>&1
        touch {output.marker}
        """
