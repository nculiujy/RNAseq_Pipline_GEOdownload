
import os
import argparse
import pandas as pd
import glob
from collections import defaultdict
import datetime

def parse_args():
    parser = argparse.ArgumentParser(description="Step 3.02: Merge Quantification Results")
    parser.add_argument("--inputdir", required=True, help="Directory containing quantification results")
    parser.add_argument("--outputdir", required=True, help="Output directory for merged matrices")
    parser.add_argument("--filter_csv", required=False, default=None, help="CSV file containing QC results (alignment_quality.csv), optional")
    # 移除 species 参数，改为自动检测
    parser.add_argument("--gtf_base", help="GTF annotation base directory")
    return parser.parse_args()

def get_species_from_path(path):
    path_lower = path.lower()
    if "homo" in path_lower:
        return "human"
    elif "mouse" in path_lower:
        return "mouse"
    return None

def get_annotation_dirs(species, gtf_base_dir):
    """Reuse annotation structure from Quant script to know which directories to look for"""
    if species == 'human':
        return [
            "mRNA/genecode/stringtie",
            "eRNA/EnhancerAtlas/stringtie",
            "eRNA/Ensembl/stringtie",
            "eRNA/FANTOM5/stringtie",
            "lncRNA/GENCODE/stringtie",
            "miRNA/miRBase/stringtie",
            "miRNA/MirGeneDB/stringtie",
            "lncRNA/NONCODE/stringtie"
        ]
    elif species == 'mouse':
        return [
            "mRNA/genecode/stringtie",
            "eRNA/EnhancerAtlas/stringtie",
            "eRNA/Ensembl/stringtie",
            "eRNA/FANTOM5/stringtie",
            "lncRNA/GENCODE/stringtie",
            "miRNA/miRBase/stringtie",
            "miRNA/MirGeneDB/stringtie",
            "lncRNA/NONCODE/stringtie"
        ]
    return []

def load_qc_pass_samples(qc_csv):
    """Load sample IDs that passed QC"""
    if qc_csv is None or not os.path.exists(qc_csv):
        print(f"Warning: QC file not provided or not found. Proceeding with all samples.")
        return None
    
    try:
        df = pd.read_csv(qc_csv)
        passed_samples = set(df[df['Passed'] == 'Yes']['Sample_ID'])
        print(f"Loaded {len(passed_samples)} passing samples from QC file.")
        return passed_samples
    except Exception as e:
        print(f"Error reading QC file: {e}")
        return None

def merge_expression_matrices(inputdir, outputdir, gtf_base, passed_samples):
    """
    Iterate through all subdirectories to find gene_abund.tab
    Detect species from path.
    """
    
    # We need to scan inputdir to find "homo" or "mouse" directories first
    # Or just search all gene_abund.tab and group by species/annotation type
    
    # Let's search for all gene_abund.tab files
    search_pattern = os.path.join(inputdir, "**", "gene_abund.tab")
    all_files = glob.glob(search_pattern, recursive=True)
    
    # Group files by (species, annotation_type)
    # File path structure: .../inputdir/homo/.../mRNA/genecode/stringtie/SampleID/gene_abund.tab
    
    grouped_files = defaultdict(list)
    
    for fpath in all_files:
        # Infer species from full path
        rel_path = os.path.relpath(fpath, inputdir)
        species = get_species_from_path(rel_path)
        
        if not species:
            continue
            
        # Infer annotation type from path
        # Check against known annotation dirs
        anno_dirs = get_annotation_dirs(species, gtf_base)
        matched_anno = None
        for anno in anno_dirs:
            if anno in rel_path:
                matched_anno = anno
                break
        
        if matched_anno:
            grouped_files[(species, matched_anno)].append(fpath)

    # Process each group
    for (species, anno_subpath), files in grouped_files.items():
        print(f"Processing {species} - {anno_subpath} ({len(files)} files)")
        
        expr_data = defaultdict(dict)
        sample_ids = set()
        
        for fpath in files:
            sample_id = os.path.basename(os.path.dirname(fpath))
            
            if passed_samples is not None and sample_id not in passed_samples:
                continue
                
            sample_ids.add(sample_id)
            
            try:
                df = pd.read_csv(fpath, sep='\t')
                for _, row in df.iterrows():
                    gene_id = row['Gene ID']
                    tpm = row['TPM']
                    expr_data[gene_id][sample_id] = tpm
            except Exception as e:
                print(f"  Error reading {fpath}: {e}")
        
        if not expr_data:
            continue
            
        matrix_df = pd.DataFrame.from_dict(expr_data, orient='index')
        matrix_df.fillna(0, inplace=True)
        matrix_df.sort_index(axis=1, inplace=True)
        matrix_df.sort_index(axis=0, inplace=True)
        
        safe_name = anno_subpath.replace('/', '_')
        # Output file: species_anno_matrix.csv
        out_file = os.path.join(outputdir, f"{species}_{safe_name}_matrix.csv")
        
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        matrix_df.to_csv(out_file)
        print(f"  Saved matrix to {out_file}")

def main():
    args = parse_args()

    # 默认使用项目内置注释目录
    if not args.gtf_base:
        args.gtf_base = "workflow/anno"

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_outputdir = args.outputdir.rstrip(os.sep)
    final_outputdir = f"{base_outputdir}_{timestamp}"
    print(f"Output directory updated with timestamp: {final_outputdir}")
    
    passed_samples = load_qc_pass_samples(args.filter_csv)
    
    merge_expression_matrices(args.inputdir, final_outputdir, args.gtf_base, passed_samples)
    
    parent_dir = os.path.dirname(args.outputdir)
    flag_file = os.path.join(parent_dir, "Merge_finished.txt")
    
    with open(flag_file, "w") as f:
        f.write(f"Merge finished at {timestamp}. Output: {final_outputdir}\n")
    print(f"Created finish flag at {flag_file}")

if __name__ == "__main__":
    main()
