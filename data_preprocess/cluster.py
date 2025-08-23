import subprocess
from pathlib import Path


def run_cd_hit(input_path: Path, output_path: Path, c=0.4, n=2, M=16000, T=8):
    cmd = [
        "cd-hit",
        "-i", str(input_path),
        "-o", str(output_path),
        "-c", str(c),
        "-n", str(n),
        "-M", str(M),
        "-T", str(T)
    ]
    subprocess.run(cmd, check=True)


def main(merged_path_str):
    merged_path = Path(merged_path_str)
    clustered_path = merged_path.with_name(merged_path.stem.replace("merged", "clustered") + merged_path.suffix)
    
    run_cd_hit(merged_path, clustered_path)

    ptm_name = merged_path.parent.name
    base_dir = merged_path.parents[1]
    no_ptm_dir = base_dir.with_name("no_ptm")
    no_ptm_path = no_ptm_dir / f"filtered_no_{ptm_name}.fasta"
    clustered_no_ptm_path = no_ptm_path.with_name(no_ptm_path.stem.replace("filtered", "clustered") + no_ptm_path.suffix)
    
    run_cd_hit(no_ptm_path, clustered_no_ptm_path)

if __name__ == "__main__":
    main()
