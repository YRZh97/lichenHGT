#!/usr/bin/env python3
import pysam
import os
import numpy as np

def mean_coverage(bam, contig, start, end, mapq=30):
    start = max(0, start)
    end = min(bam.get_reference_length(contig), end)

    if start >= end:
        return 0
    # count_coverage 默认 0-based, 左闭右开
    cov = bam.count_coverage(contig, start, end, quality_threshold=mapq)
    # cov 是 (A,C,G,T) 四个数组，长度 = end-start
    total_cov = np.sum(cov, axis=0)  # 每个位点总覆盖度
    return np.mean(total_cov) if len(total_cov) > 0 else 0

def depth(bamfile,contig,pos,window,gene):
    sample = os.path.basename(bamfile).replace(".rmdup.bam", "")
    bam = pysam.AlignmentFile(bamfile, "rb")

    depth = mean_coverage(bam, contig, pos-1, pos)  # pysam 用 0-based 左闭右开
    left_coverage = mean_coverage(bam, contig, pos-window, pos-1)
    right_coverage = mean_coverage(bam, contig, pos, pos+window)

    start = max(0,pos - window)
    end = min(pos + window,bam.get_reference_length(contig)) 
    crossing_reads = 0

    for read in bam.fetch(contig, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < 30:
            continue
        read_start = read.reference_start
        read_end = read.reference_end
        if read_start <= (pos-1) and read_end >= pos:
            if read.has_tag('NM') and read.get_tag('NM') == 0:
                crossing_reads += 1
    bam.close()
    return [sample, gene, depth, left_coverage, right_coverage, crossing_reads]

def getpos(mag,gene):
    pos = []
    contig = None
    gff_file = 'gff/' + mag + '.gff'
    if not os.path.exists(gff_file):
        return None, None
    for line in open(gff_file,'r'):
        if '#' in line:
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        if fields[2] == 'gene' and gene in fields[-1]:
            contig = fields[0]
            pos = [int(fields[3]), int(fields[4])]
            break
    return contig,pos if contig else (None, None)

if __name__ == "__main__":
    output = open('depth.txt','w')
    output.write("Sample\tGene\tContig\tBreakpoint_cov\tLeft_cov\tRight_cov\tCrossing_reads\tOrder\tPos\n")
    for line in open('allHGTs.txt','r'):#allHGTs.txt: a namelist of all HGT genes
        genome = line.strip().rsplit('_',1)[0]
        gene = line.strip().rsplit('_',1)[1][:-3]
        sampleid = genome.split('_')[0]
        contig,pos = getpos(genome,gene)
        output.write(f"{sampleid}\t{gene}\t{contig}\t")

        if not contig or not pos:
            output.write(f"No GFF match\n")
            continue
        
        bamfile = f"breakpoint/{sampleid}/{sampleid}.rmdup.bam"
        if os.path.exists(bamfile):
            for i in range(len(pos)):
                #print(f"{sampleid},{gene},{p}")
                result = depth(bamfile,contig,pos[i],150,gene)
                newline = '\t'.join(map(str, result)) + f"\t{str(i)}\t{str(pos[i])}\n"
                output.write(newline)
    output.close()


