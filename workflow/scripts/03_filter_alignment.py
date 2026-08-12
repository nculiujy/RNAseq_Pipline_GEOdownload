
import os
import argparse
import csv
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Step 3.3: Filter Alignment Quality")
    parser.add_argument("--inputdir", required=True, help="Directory containing results")
    parser.add_argument("--outputdir", required=True, help="Output directory for reports")
    parser.add_argument("--cutoff", type=float, default=70.0,
                        help="Alignment rate cutoff (default: 70.0%%)")
    return parser.parse_args()

def extract_alignment_rates(inputdir, outputdir, cutoff=70.0):
    print(f"Extracting alignment rates (cutoff={cutoff}%)...")
    results = []
    
    # 递归查找所有 QC_results.log
    for root, dirs, files in os.walk(inputdir):
        if "QC_results.log" in files:
            log_file = os.path.join(root, "QC_results.log")
            # 假设路径结构 .../GSExxx/hisat2file/SampleID/QC_results.log
            # 向上两级是 GSE 目录，父目录是 SampleID
            sample_id = os.path.basename(root)
            
            try:
                with open(log_file, "r") as f:
                    content = f.read()
                    matches = re.findall(r"(\d+\.\d+)%", content)
                    if matches:
                        rate = float(matches[-1])
                        passed = "Yes" if rate >= cutoff else "No"
                        results.append({
                            "Sample_ID": sample_id,
                            "Alignment_Rate": rate,
                            "Passed": passed,
                            "Path": log_file
                        })
            except Exception as e:
                print(f"Failed to read log {log_file}: {e}")
    
    csv_file = os.path.join(outputdir, "alignment_quality.csv")
    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Sample_ID", "Alignment_Rate", "Passed", "Path"])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Written {len(results)} records to {csv_file} (cutoff={cutoff}%)")
    return csv_file

def main():
    args = parse_args()
    os.makedirs(args.outputdir, exist_ok=True)
    csv_file = extract_alignment_rates(args.inputdir, args.outputdir, cutoff=args.cutoff)
    
    with open(os.path.join(args.outputdir, "Filter_finished.txt"), "w") as f:
        f.write(f"Filtering finished (cutoff={args.cutoff}%).\n")

if __name__ == "__main__":
    main()
