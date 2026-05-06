import os
from collections import defaultdict
import re
from Bio import SeqIO

# Configuration
NAME_FILE = "hq.name"  # MAG\taxonomy information
OUT_FILE = "merged_top_hit.tsv"
GFF_DIR = "gfffile"  # Directory containing all GFF files
FASTA_FILE = "hgt.fa"  # Sequences of candidate HGTs
FASTA_DIRS = "path/to/genomefiles"  # Directory with genome FASTA files
TAXID_FILE = "taxid"  # TaxIDs for target group
CONTIG_OUT_FILE = "gene_on_contig"

def normalize_id(name):
    """Normalize ID by removing suffixes."""
    name = name.replace(".gff", "").replace(".fasta", "").replace(".fna", "")
    name = re.sub(r"\.fa_bin$", "", name)
    return name

def get_genome_id(name):
    """Extract genome ID from name."""
    m = re.match(r"(GCA_\d+\.\d+)", name)
    if m:
        return m.group(1)
    m = re.match(r"(.+?_Refined_\d+)", name)
    if m:
        return m.group(1)
    return None

def extract_top_hits():
    """Extract top hits from Diamond output files."""
    with open(NAME_FILE) as f:
        lines = [line.strip() for line in f if line.strip()]

    with open(OUT_FILE, "w") as out:
        for line in lines:
            parts = line.split(",")
            faa_name = parts[0]
            sample_dir = faa_name.split('_')[0]
            prefix = faa_name.replace(".fa_bin.faa", "")
            diamond_file = os.path.join(sample_dir, "max300_subseq.fmt6")

            if not os.path.exists(diamond_file):
                print(f"Missing: {diamond_file}")
                continue

            seen = set()
            print(f'Processing: {faa_name}')
            with open(diamond_file) as f:
                for line in f:
                    cols = line.strip().split("\t")
                    qseqid = cols[0]
                    sseqid = cols[1]
                    taxid = cols[-1]

                    if qseqid in seen:
                        continue
                    seen.add(qseqid)

                    new_query = f"{faa_name}_{qseqid[:-3]}"
                    out.write(f"{new_query}\t{sseqid}\t{taxid}\n")

def load_gene2taxid():
    """Load gene to taxid mapping."""
    gene2taxid = {}
    with open(OUT_FILE) as f:
        next(f)  # Skip header if present
        for line in f:
            cols = line.strip().split("\t")
            query = cols[0]
            taxid = cols[-1].split(';')[0]
            gene2taxid[query] = taxid
    return gene2taxid

def load_target_taxids():
    """Load target taxids."""
    target_taxids = set()
    with open(TAXID_FILE) as f:
        for line in f:
            target_taxids.add(line.strip())
    return target_taxids

def collect_target_samples():
    """Collect target samples and contig lengths."""
    target_samples = set()
    contig_length = defaultdict(int)

    for f in os.listdir(FASTA_DIRS):
        if f.endswith((".fasta", ".fna")):
            sample_id = normalize_id(f)
            target_samples.add(sample_id)
            genome_id = get_genome_id(sample_id)
            print(f"Processing: {f}...")
            for record in SeqIO.parse(os.path.join(FASTA_DIRS, f), "fasta"):
                contig_full = f"{genome_id}|{record.id}"
                contig_length[contig_full] = len(record.seq)
    print(f"=====Total target samples: {len(target_samples)}=====")
    return target_samples, contig_length

def collect_gff_info(target_samples):
    """Collect GFF information."""
    gene2contig = {}
    contig_gene_count = defaultdict(int)
    contig2genes = defaultdict(set)

    for file in os.listdir(GFF_DIR):
        if not file.endswith(".gff"):
            continue
        print(f'Processing: {file}')
        gff_path = os.path.join(GFF_DIR, file)
        genome_id = normalize_id(file)
        if genome_id not in target_samples:
            continue

        genome_id = get_genome_id(genome_id)
        with open(gff_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 9 or cols[2] != "gene":
                    continue

                contig_full = f"{genome_id}|{cols[0]}"
                contig_gene_count[contig_full] += 1
                gene_id = cols[8].strip()
                key = f"{genome_id}_{gene_id}"
                gene2contig[key] = contig_full
                contig2genes[contig_full].add(key)

    return gene2contig, contig_gene_count, contig2genes

def parse_fasta_target_genes(gene2contig):
    """Parse FASTA file for target genes."""
    contig_target_genes = defaultdict(set)
    unmatched_genes = []

    with open(FASTA_FILE) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line.strip().replace(">", "").replace("eu_", "")
            gene_match = re.search(r"(g\d+)", header)
            genome_match = re.search(r"^[^_]+_(.+?)_g\d+", header)
            if not gene_match or not genome_match:
                continue
            genome_id = genome_match.group(1)
            key = f"{genome_id}_{gene_match.group(1)}"
            contig = gene2contig.get(key)
            if contig:
                contig_target_genes[contig].add(key)
            else:
                unmatched_genes.append(key)

    return contig_target_genes, unmatched_genes

def write_contig_stats(contig_length, contig_gene_count, contig_target_genes, contig2genes, gene2taxid, target_taxids):
    """Write contig statistics to output file."""
    with open(CONTIG_OUT_FILE, "w") as out:
        out.write("Contig\tTotal_gene\tTarget_gene\tRatio\tTarget_gene_list\t20k\tnonHGT_taxids\tnonHGT_fungi\tnonHGT_nonfungi\n")
        for contig in contig_length:
            total = contig_gene_count.get(contig, 0)
            target = len(contig_target_genes.get(contig, set()))
            ratio = target / total if total > 0 else 0
            gene_list = ",".join(sorted(contig_target_genes.get(contig, set())))
            judge = 1 if contig_length[contig] > 20000 else 0

            taxid_str = ""
            target_count = ""
            non_target_taxids = ""
            hgt_genes = contig_target_genes.get(contig, set())

            if hgt_genes:
                all_genes = contig2genes.get(contig, set())
                non_hgt_genes = all_genes - hgt_genes
                taxids = []
                nonfungi = []
                target_count_val = 0

                for g in non_hgt_genes:
                    taxid = gene2taxid.get(g)
                    if not taxid:
                        nonfungi.append("NA")
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

def main():
    # Step 1: Extract top hits
    extract_top_hits()

    # Step 2: Load mappings
    gene2taxid = load_gene2taxid()
    target_taxids = load_target_taxids()

    # Step 3: Collect samples and contigs
    target_samples, contig_length = collect_target_samples()

    # Step 4: Collect GFF info
    gene2contig, contig_gene_count, contig2genes = collect_gff_info(target_samples)

    # Step 5: Parse target genes
    contig_target_genes, unmatched_genes = parse_fasta_target_genes(gene2contig)

    # Step 6: Write statistics
    write_contig_stats(contig_length, contig_gene_count, contig_target_genes, contig2genes, gene2taxid, target_taxids)

    print("Pipeline completed successfully!")

if __name__ == "__main__":
    main()
