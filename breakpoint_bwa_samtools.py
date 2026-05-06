#!/usr/bin/env python3

import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

#####CAUTION#####
# conda run -n {env} may not work in some environments, you might need to activate the env first and run the commands directly.

bwa_env = "bwa"
samtools_env = "samtools"

threads = 10
parallel_jobs = 2

def run_pipeline(genome):
    prefix = genome.split('_')[0]
    sample_dir = os.path.join('breakpoint',prefix)
    os.makedirs(sample_dir, exist_ok=True)
    print(f"=== Processing {prefix} ===")

    cwd = os.getcwd()
    os.chdir(sample_dir)
    try:
    # 自动寻找对应的 reads
        reads_R1 = f"sample_fq/{prefix}_1.fastq.gz"
        reads_R2 = f"sample_fq/{prefix}_2.fastq.gz"
        if not (os.path.exists(reads_R1) and os.path.exists(reads_R2)):
            return f"{prefix} failed: {reads_R1} or {reads_R2} not found"

        cmd_index = f"conda run -n {bwa_env} bwa index samples/{prefix}/final.contigs.fa -p {prefix}_contigs"
        #subprocess.run(cmd_index, shell=True, check=True)
        
        cmd_mem =  f"bwa mem -t {threads} {prefix}_contigs {reads_R1} {reads_R2} > {prefix}.sam"
        subprocess.run(cmd_mem, shell=True, check=True)

        cmd_samtools_view = f"conda run -n {samtools_env} samtools view -bS {prefix}.sam > {prefix}.bam"
        #subprocess.run(cmd_samtools_view, shell=True, check=True)

        name_sorted = f"{prefix}_sorted.bam"
        cmd_sort_name = f"conda run -n {samtools_env} samtools sort -n -@ {threads} -o {name_sorted} {prefix}.bam"
        #subprocess.run(cmd_sort_name, shell=True, check=True)

        fixmate_bam = f"{prefix}.fixmate.bam"
        cmd_fixmate = f"conda run -n {samtools_env} samtools fixmate -m {name_sorted} {fixmate_bam}"
        #subprocess.run(cmd_fixmate, shell=True, check=True)

        pos_sorted = f"{prefix}.pos_sorted.bam"
        cmd_sort_pos = f"conda run -n {samtools_env} samtools sort -@ {threads} -o {pos_sorted} {fixmate_bam}"
        #subprocess.run(cmd_sort_pos, shell=True, check=True)

        rmdup_bam = f"{prefix}.rmdup.bam"
        cmd_markdup = f"conda run -n {samtools_env} samtools markdup -r -@ {threads} {pos_sorted} {rmdup_bam}"
        #subprocess.run(cmd_markdup, shell=True, check=True)
        
        cmd_index_bam = f"conda run -n {samtools_env} samtools index {rmdup_bam}"
        #subprocess.run(cmd_index_bam, shell=True, check=True)

        return f"{prefix} finished, result: {rmdup_bam}"

    except subprocess.CalledProcessError as e:
        return f"{prefix} failed: {e}"

    finally:
        os.chdir(cwd)

if __name__ == "__main__":
    genomes = []
    for line in open('hq.name','r'):
        genomes.append(line.strip().split('.')[0])

    results = []
    with ProcessPoolExecutor(max_workers=parallel_jobs) as executor:
        future_to_genome = {executor.submit(run_pipeline, genome): genome for genome in genomes}
        for future in as_completed(future_to_genome):
            res = future.result()
            print(res)


