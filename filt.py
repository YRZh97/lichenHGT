# !/user/bin/env python
# -*- coding:utf-8 -*-

import os,re,shutil
from Bio import SeqIO

def makedirectory(path,filename):
    absolute_path = path + '/' + filename
    if os.path.exists(absolute_path):
        shutil.rmtree(absolute_path)  # 删除目录，包括目录下的所有文件
    os.makedirs(absolute_path)
    return
    
'''
#this step has been abandoned.
badcontig = set()
with open("kmerfreq_total_4mer") as f:#kmerfreq: samplename\tcontig1,contig2
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) > 1:
            prefix = parts[0]
            for k in parts[1].split(","):
                badcontig.add(f"{prefix}|{k}")
'''

names = []#hq.name: MAG\ttaxonomy information
for line in open('hq.name','r'):
    names.append(line.strip().split(',')[0])

for n in names:
    if os.path.exists(f'newdir/{n}'):
        pass
    else:
        makedirectory('newdir', n)
        print(f'{f} dir successfully created')
        
        t = 0
        p = []#contigs longer than 20k and not excluded based on 4-mer frequency
        for record in SeqIO.parse(f'bin/{n}_genomic.fna', 'fasta'):
            t += 1
            if len(record.seq) <= 20000 or record.id in badcontig:
                pass
            else:
                p.append(record.id)
        print(t,len(p))

        gffdic = {}
        for line in open(f'gff/{n}.gff', 'r'):
            if '\ttranscript\t' in line:
                gffdic[line.strip().split('\t')[-1]] = line.strip().split('\t')[0]#AUGUSTUS gff
        print(len(gffdic))

        m = 0
        outputfile = open(f'newdir/{n}.filtered.faa','w')
        for rec in SeqIO.parse(f'faa/{n}_protein.faa', 'fasta'):
            if gffdic[rec.id] in p:
                SeqIO.write(rec, outputfile, 'fasta')
                m += 1
            else:
                pass
        print(m)
