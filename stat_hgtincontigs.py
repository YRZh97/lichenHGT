from collections import defaultdict
import re
import os
from Bio import SeqIO

gff_dir = "gfffile"        # gff directory containing all gff files
fasta_file = "hgt.fa"   # sequences of candidate HGTs


fasta_dirs = 'path/to/genomefiles'#

def normalize_id(name):
    # 去掉后缀
    name = name.replace(".gff", "").replace(".fasta", "").replace(".fna", "")
    
    # 🟢 去掉 Refined 类型的 .fa_bin
    name = re.sub(r"\.fa_bin$", "", name)
    
    return name

def get_genome_id(name):
    # GCA类型
    m = re.match(r"(GCA_\d+\.\d+)", name)
    if m:
        return m.group(1)

    # 数字Refined类型
    m = re.match(r"(.+?_Refined_\d+)", name)
    if m:
        return m.group(1)

    return None


target_samples = set()
contig_length = defaultdict(int)

gene2taxid = {}   # key: genome_gene  → taxid
target_taxids = set()  # 你的目标类群taxid集合

diamond_file = 'merged_top_hit.tsv'

with open(diamond_file) as f:
    next(f)  # 跳过表头
    for line in f:
        cols = line.strip().split("\t")
        
        query = cols[0]   # 20211966_Refined_1|g1.t1
        taxid = cols[-1].split(';')
        gene2taxid[query] = taxid[0]
#print(gene2taxid)

with open("taxid") as f:#taxids for target group, extracted with 'taxonkit list' command
    for line in f:
        target_taxids.add(line.strip())
#print(target_taxids)


for f in os.listdir(fasta_dirs):
    if f.endswith((".fasta", ".fna")):
    #if f == 'GCA_964257355.1_private_T1894_concoct_bin.66_genomic.fna':
        sample_id = normalize_id(f)

        target_samples.add(sample_id)
        genome_id = get_genome_id(sample_id)
        print(f"processing:{f}...")
        for record in SeqIO.parse(f'{d}/{f}', "fasta"):
            contig = record.id
            contig_full = f"{genome_id}|{contig}"
            contig_length[contig_full] = len(record.seq)
print(f"=====Total target samples: {len(target_samples)}=====")

gene2contig = {}
contig_gene_count = defaultdict(int)
contig2genes = defaultdict(set)

for file in os.listdir(gff_dir):
    if not file.endswith(".gff"):
        continue
    print(f'processing:{file}')
    gff_path = os.path.join(gff_dir, file)
    genome_id = normalize_id(file)
    if genome_id not in target_samples:
        continue
    
    genome_id = get_genome_id(genome_id)
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            
            contig = cols[0]
            feature = cols[2]
            attr = cols[8]
            
            if feature == "gene":
                contig_full = f"{genome_id}|{contig}"  # 防止不同gff contig重名
                contig_gene_count[contig_full] += 1
                
                gene_id = attr.strip()  # g973
                
                key = f"{genome_id}_{gene_id}"
                #print(f'gff-based key:{key}')
                gene2contig[key] = contig_full
                contig2genes[contig_full].add(key)


contig_target_genes = defaultdict(set)
unmatched_genes = []

with open(fasta_file) as f:
    for line in f:
        if not line.startswith(">"):
            continue
        
        header = line.strip().replace(">", "")
        header = header.replace("eu_", "")
        # 🔥 通用解析（最关键）
        # 匹配最后的 gXXX
        gene_match = re.search(r"(g\d+)", header)
        if not gene_match:
            continue
        gene_id = gene_match.group(1)
        
        # genome_id = 去掉前缀和 _gXXX.t1
        genome_match = re.search(r"^[^_]+_(.+?)_g\d+", header) 
        if not genome_match: 
            continue 
        genome_id = genome_match.group(1)
        
        key = f"{genome_id}_{gene_id}"
        #print(f'fasta-based key:{key}')
        if key in gene2contig:
            contig = gene2contig[key]
            contig_target_genes[contig].add(key)
        else:
            unmatched_genes.append(key)

out_file = "gene_on_contig"
with open(out_file, "w") as out:
    out.write("Contig\tTotal_gene\tTarget_gene\tRatio\tTarget_gene_list\t20k\tnonHGT_taxids\tnonHGT_fungi\tnonHGT_nonfungi\n")
#print(f"\nUnmatched genes: {len(unmatched_genes)}")
#for g in unmatched_genes[:50]:   # 只看前50个，避免刷屏
    #print(g)
    for contig in contig_length:
        total = contig_gene_count.get(contig, 0)
        target = len(contig_target_genes.get(contig,set()))
        ratio = target / total if total > 0 else 0 
        gene_list = ",".join(sorted(contig_target_genes.get(contig,set())))
        judge = 1 if contig_length[contig] > 20000 else 0

        taxid_str = ""
        target_count = ""
        non_target_taxids = ""
        hgt_genes = contig_target_genes.get(contig,set())

        if len(hgt_genes) > 0:
            all_genes = contig2genes.get(contig,set())
            non_hgt_genes = all_genes - hgt_genes
            taxids = []
            nonfungi = []
            target_count_val = 0

            for g in non_hgt_genes:
                taxid = gene2taxid.get(g)

                if not taxid:
                    nonfungi.append("NA")   # 或者 "NoHit"
                    continue

                taxids.append(taxid)

                if taxid in target_taxids:
                    target_count_val += 1
                else:
                    nonfungi.append(taxid)

            taxid_str = ";".join(taxids)
            target_count = str(target_count_val)
            non_target_taxids = ";".join(nonfungi)

        out.write(f"{contig}\t{total}\t{target}\t{ratio:.4f}\t{gene_list}\t{judge}\t{taxid_str}\t{target_count}\t{non_target_taxids}\n")
