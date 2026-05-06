#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os

BASE_PATH = "path/to/directory" 
NR_DB = "nr_diamond_taxid.dmnd"#diamond database for NR with taxid information, using --taxonmap, --taxonnodes and --taxonnames parameters for 'diamond makedb' command
EU_DB = "90_5.dmnd"#diamond database for high-completeness MAGs generated in this study with default parameters for 'diamond makedb' command

# 读取样本列表
with open("hq.name") as f:
    samples = [line.strip().split(",")[0] for line in f]

for sample in samples:
    sample_prefix = sample.split("_")[0]
    sample_dir = os.path.join(BASE_PATH, sample_prefix)
    files = os.listdir(sample_dir)
    
    # 只处理目录中只有一个文件的情况
    if len(files) != 1:
        continue
    
    # 查找 filtered.faa 文件
    for filename in files:
        if not filename.endswith("filtered.faa"):
            continue
        
        print(sample)
        print(filename)
        
        faa_path = os.path.join(sample_dir, filename)
        
        # 命令1: NR 数据库搜索
        diamond1 = (
            f"diamond blastp --db {NR_DB} -q {faa_path} "
            f"-o {os.path.join(sample_dir, 'max300_subseq.fmt6')} "
            f"--max-target-seqs 300 --outfmt 6 qseqid sseqid pident length mismatch "
            f"gapopen qstart qend sstart send evalue bitscore full_sseq staxids --threads 5"
        )
        print(diamond1)
        # os.system(diamond1)
        
        # 命令2: 自定义数据库搜索
        output_name = f"eu_{filename.split('.')[0]}.fmt6"
        diamond2 = (
            f"diamond blastp --db {EU_DB} -q {faa_path} "
            f"-o {os.path.join(sample_dir, output_name)} "
            f"--outfmt 6 --threads 5"
        )
        print(diamond2)
        # os.system(diamond2)

