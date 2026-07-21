#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
from collections import defaultdict
from ete3 import Tree

BASE_DIR = "/path/to/your/directory"  
TREEFILE = "treefile_hgt"#name list of HGTs, format: sample_id\tgene_id
HQ_FILE = "hq.name"
FINAL_OUTPUT = os.path.join(BASE_DIR, "donor_species_final")
TREES_DIR = os.path.join(BASE_DIR, "treefile")
HQMAG_TAXON_FILE = os.path.join(BASE_DIR, "hqmag.taxon")#taxonomy information of high-completeness MAGs, format: genome_id\tlineage


def lca(lineage1, lineage2):
    common = []
    for a, b in zip(lineage1.split(";"), lineage2.split(";")):
        if a == b and a.strip():
            common.append(a)
        else:
            break
    return ";".join(common) if common else "donorconfusion"


def read_tsv(path, sep="\t", min_cols=0):
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split(sep)
            if len(parts) >= min_cols:
                yield parts


def load_dict(path, key_index=0, value_index=1, sep="\t", min_cols=2):
    return {cols[key_index]: cols[value_index] for cols in read_tsv(path, sep=sep, min_cols=min_cols)}


def build_namedic(path):
    namedic = defaultdict(list)
    with open(path) as f:
        for line in f:
            sample_id, gene_id = line.split()
            namedic[sample_id].append(gene_id)
    return namedic


def find_top_candidate(gene, fungilis):
    path = os.path.join(TREES_DIR, f"{gene}_new_tree.nw")
    try:
        tree = Tree(path, format=2, quoted_node_names=True)
    except Exception:
        tree = Tree(path, format=5, quoted_node_names=True)

    node = tree.search_nodes(name=gene)[0]
    while node:
        if node.is_leaf():
            node = node.up
            continue
        for leaf in node:
            name = leaf.name
            if name.split(".")[0] in fungilis or "Viridiplantae" in name:#“Viridiplantae” is used for proteins in algal MAGs; for proteins in fungal MAGs, “Fungi” is used to determine their taxonomic origin.
                continue
            subject = name.split("|")[0]
            return subject
        node = node.up
    return None


def update_identities(path, identdic, taxon=None):
    for cols in read_tsv(path, min_cols=3):
        query, subject, ident = cols[0], cols[1], float(cols[2])
        if query in identdic and subject in identdic[query]:
            identdic[query][subject] = ident
            if taxon is not None:
                taxon[subject] = cols[-1]

namedic = build_namedic(TREEFILE)
fungilis = [line.split(",")[1] for line in open(HQ_FILE)]
taxon3 = load_dict(HQMAG_TAXON_FILE, 0, 1)

with open(FINAL_OUTPUT, "w") as output2:
    output2.write("#sample\tgene\tMMSH_iden\tMMSH\tMMSH_taxid\tdonor\tbbhO\n")

    for sample_id, genes in namedic.items():
        taxon2 = load_dict(os.path.join(BASE_DIR, sample_id, "reformat.txt"), 0, 1)#reformat.txt contains taxonomic information of matched proteins, format: subject_id\ttaxonomy lineage
        donor_path = os.path.join(BASE_DIR, sample_id, "donor_species")

        newdic = {}
        identdic = {}
        for gene in genes:
            newdic[gene] = []
            identdic[gene] = {}
            print(gene)
            subject = find_top_candidate(gene, fungilis)
            if subject:
                newdic[gene].append(subject)
                identdic[gene][subject] = 0.0

        update_identities(os.path.join(BASE_DIR, sample_id, f"{sample_id}.fmt6"), identdic)#sample_id.fmt6 contains the alignment results against the newly assembled MAGs obtained from sequencing data generated in this study.
        taxon_map = {}
        update_identities(os.path.join(BASE_DIR, sample_id, "max300_subseq.fmt6"), identdic, taxon_map)#max300_subseq.fmt6 contains the alignment results against NR database

        with open(donor_path, "w") as totaloutput:
            for gene, hits in identdic.items():
                best_hit, best_score = "error", 0
                for subject, score in hits.items():
                    if score > best_score:
                        best_hit, best_score = subject, score

                if best_hit.startswith("GCA_"):
                    taxid = taxon3.get(best_hit.split(".")[0])
                    if taxid:
                        totaloutput.write(
                            f"{sample_id}\t{gene}\t{best_score}\t{best_hit}\t-\t{taxid}\n"
                        )
                    else:
                        totaloutput.write(f"{sample_id}\t{gene}\t{best_score}\t{best_hit}\n")
                else:
                    taxon_value = taxon_map[best_hit]
                    genus_id = taxon_value.split(";")[0]
                    if genus_id in taxon2:
                        totaloutput.write(
                            f"{sample_id}\t{gene}\t{best_score}\t{best_hit}\t{taxon_value}\t{taxon2[genus_id]}\n"
                        )
                    else:
                        totaloutput.write(f"{sample_id}\t{gene}\t{best_score}\t{best_hit}\t{taxon_value}\n")

        revised_path = os.path.join(BASE_DIR, sample_id, "donor_species_revised")
        bbho = {}
        for cols in read_tsv(os.path.join(BASE_DIR, sample_id, f"{sample_id}_ai.txt")):#{sample_id}_ai.txt was generated in the ai_calc.py script
            query = cols[0]
            if query not in namedic[sample_id]:
                continue
            subject = cols[5]
            if subject.startswith("GCA_"):
                bbho[query] = [";".join(taxon3[subject.split(".")[0]].split(";")[:3]), subject]
            elif subject.split("|")[1] in taxon2:
                bbho[query] = [";".join(taxon2[subject.split("|")[1]].split(";")[:3]), subject]

        with open(os.path.join(BASE_DIR, sample_id, "donor_species")) as infile, open(revised_path, "w") as output:
            output.write("#sample\tgene\tMMSH_iden\tMMSH\tMMSH_taxid\tMMSH_donor\tbbhO_donor\tbbhO\n")
            for cols in read_tsv(os.path.join(BASE_DIR, sample_id, "donor_species")):
                if len(cols) < 7:
                    cols += [""] * (7 - len(cols))
                if cols[0] != sample_id or cols[1] not in bbho:
                    continue
                donor_name, subject = bbho[cols[1]]
                relation = "consis" if donor_name in cols[-1] else "noconsis"
                output.write("\t".join(cols) + f"\t{donor_name}\t{subject}\t{relation}\n")

        for cols in read_tsv(revised_path):
            status = cols[-1]
            if status == "consis":
                output2.write("\t".join(cols[:-3] + [cols[-2]]) + "\n")
            elif status == "noconsis":
                if cols[-3] in cols[5]:
                    output2.write("\t".join(cols[:-3] + [cols[-2]]) + "\n")
                else:
                    output2.write(
                        "\t".join(cols[:5] + [lca(cols[5], cols[6]), cols[-2]]) + "\n"
                    )

