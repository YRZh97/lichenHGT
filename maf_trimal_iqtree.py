#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import multiprocessing

def build_commands(filename):
    if filename.endswith('_prealign.fa'):
        base = filename.rsplit('_', 1)[0]
        mafft_cmd = f'mafft --thread 3 --auto {filename} > {base}_prealign.maf'
        trimal_cmd = f'trimal -in {base}_prealign.maf -out {base}_prealign.trimal -matrix newmatrix -keepheader -automated1'
        iqtree_cmd = f'iqtree --redo -nt 5 -st AA -s {base}_prealign.trimal -m TEST -mrate G4 -keep-ident -bb 1000 -pre {base}'
        return {
            'mafft': mafft_cmd,
            'trimal': trimal_cmd,
            'iqtree': iqtree_cmd,
            'base': base
        }
    return None

def run_mafft(filename):

    cmd_dict = build_commands(filename)
    if cmd_dict and os.path.exists(filename):
        print(cmd_dict['mafft'])
        os.system(cmd_dict['mafft'])

def run_trimal(filename):

    cmd_dict = build_commands(filename)
    if cmd_dict and os.path.exists(cmd_dict['base'] + '_prealign.maf'):
        print(cmd_dict['trimal'])
        os.system(cmd_dict['trimal'])

def run_iqtree(filename):

    cmd_dict = build_commands(filename)
    if cmd_dict and os.path.exists(cmd_dict['base'] + '_prealign.trimal'):
        print(cmd_dict['iqtree'])
        os.system(cmd_dict['iqtree'])

def run_pipeline(filename):

    cmd_dict = build_commands(filename)
    if not cmd_dict:
        return
    
    # 运行 MAFFT
    if os.path.exists(filename):
        print(cmd_dict['mafft'])
        os.system(cmd_dict['mafft'])
    
    # 运行 TriML
    if os.path.exists(cmd_dict['base'] + '_prealign.maf'):
        print(cmd_dict['trimal'])
        os.system(cmd_dict['trimal'])
    
    # 运行 IQ-TREE
    if os.path.exists(cmd_dict['base'] + '_prealign.trimal'):
        print(cmd_dict['iqtree'])
        os.system(cmd_dict['iqtree'])

def get_sequence_files(fpath):
    """获取所有 _prealign.fa 文件列表"""
    files = []
    for f in os.listdir(fpath):
        if f.endswith('_prealign.fa'):
            files.append(f)
    return files

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <path> [step]")
        print("  step: mafft | trimal | iqtree | all (default: all)")
        sys.exit(1)
    
    fpath = sys.argv[1]
    step = sys.argv[2].lower() if len(sys.argv) > 2 else 'all'
    
    os.chdir(fpath)
    
    files = get_sequence_files(fpath)
    print(f"Found {len(files)} sequence files: {files}")
    
    # 选择处理函数
    if step == 'mafft':
        process_func = run_mafft
    elif step == 'trimal':
        process_func = run_trimal
    elif step == 'iqtree':
        process_func = run_iqtree
    else:  # 'all' 或默认
        process_func = run_pipeline
    
    # 多进程处理
    pool_size = 3 if step in ['mafft', 'trimal'] else 5 if step == 'iqtree' else 3
    p = multiprocessing.Pool(pool_size)
    p.map(process_func, files)
    p.close()
    p.join()
    
    print("All tasks completed!")

