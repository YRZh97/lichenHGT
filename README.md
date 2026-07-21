# Codes for HGT detection in lichen metagenomes

map.R -------------------- corresponding to Fig.1  
process_fcs_res.py ------- flag contaminant sequences, corresponding to 'Contig filtering and contamination detection' part in Fig.2  
filt.py ------------------ select sequences based on the results of the previous two scripts, corresponding to 'Contig filtering and contamination detection' part in Fig.2  
diamond.py --------------- sequence alignment, corresponding to 'Alien Index score calculation' part in Fig.2  
ai_calc.py --------------- AI calculation, corresponding to 'Alien Index score calculation' part in Fig.2  
maf_trimal_iqtree.py ----- tree construction based on AI screening, corresponding to 'Tree inference' part in Fig.2  
stat_hgtincontigs.py ----- assign confidence levels to candidate genes, corresponding to 'Confidence classification' part in Fig.2  
donor_species.py --------- donor inference, corresponding to 'Donor inference' part in Fig.2  
breakpoint_bwa_samtools.py obtain coverage information at the breakpoint, corresponding to Fig. 5b &5c  
depth.py ----------------- coverage calculation, corresponding to Fig. 5b&5c

=========================================================================================================================================================================================================================  
lichenHGT-plots.ipynb contains codes used to generate the figures in the main text and supplementary materials. 'breakp_data_update.xlsx' file could be generated using depth.py, 'target_gene_annotation_2021.tsv' file is a combination of 'gene_on_contig0422' and Supplementary Datasets 4-6, 'supplementary tables1201-github.xlsx' contains fungal HGT-acquired genes and eggNOG-mapper results for total 2,021 HGT candidates.




## Step 1. Contig filtering and contamination detection
This step filters contigs by length and removes potentially contaminated contigs and genes. The primary scripts used are `filt.py` and `process_fcs_res.py`.

### Input Preparation

Place the input files into separate directories according to their file types. For example:

- Genome sequences in the `bin` directory.
- GFF annotation files in the `gff` directory.
- CDS sequence files in the `cds` directory.
- Protein sequence files in the `faa` directory.

Prepare a file named `hq.name` containing the following information for each MAG:

- MAG file name
- Genome accession
- Scientific name
- NCBI Taxonomy ID

For example:

```text
public_SRR14721965_metabat2_bin.1,GCA_964256395,Cladonia petrophila,195777
```

The taxonomy ID can be obtained from the scientific name using the `name2taxid` command in TaxonKit.

Prepare a file named `treefile` containing genome accession and geneID for each gene.

For example:

```test
GCA_964256395,g353
```

### Filtering

Run `FCS-GX` on both the genome sequence and the corresponding CDS sequence for each MAG to identify potential contaminant sequences.

For example:

```text
python3 fcs.py screen genome --fasta bin/{name}_genomic.fna --out-dir fcsgx_test/{name}_genome --gx-db '/path/to/gxdb' --tax-id {taxid}
```

After running `FCS-GX` for all MAGs, merge all generated `{n}.fcs_gx_report.txt` files into a single file named `totalfcsgx.txt`.

Like this:

```text
#genome level   seq_id  start_pos       end_pos seq_len action  div     agg_cont_cov    top_tax_name
GCA_964257225   contig  CAXVDQ010000002.1       1       1808    590398  TRIM    virs:prokaryotic viruses        55      Heterodera schachtii
```

> **Note:** Regardless of the action in the fifth column (*The recommended action*) of the `{n}.fcs_gx_report.txt`, all sequences flagged by `FCS-GX` are treated uniformly.

This merged report serves as the input for `process_fcs_res.py` in the subsequent contamination filtering step.

Then run:

```bash
python3 filt.py
```

The script will create a new directory named `newdir` in the current directory. The filtered protein sequences for each MAG will be saved as:

```text
newdir/{n}.filtered.faa
```

where `{n}` represents the corresponding MAG name.

```bash
python3 process_fcs_res.py
```

The sequences reported by the script are considered potential contaminants and are excluded from all subsequent analyses.

## Step 2. Alien Index Score Calculation

### Build the DIAMOND Database

When building the DIAMOND database, include the taxonomy-related options `--taxonmap`, `--taxonnodes`, and `--taxonnames` to enable taxonomic annotation of alignment results.

```bash
diamond makedb \
    --in nr.fa \
    --db database_name \
    --taxonmap prot.accession2taxid.gz \
    --taxonnodes nodes.dmp \
    --taxonnames names.dmp
```

### Align sequences using DIAMOND

```bash
python3 diamond.py
```

This script uses the filtered protein sequences (`newdir/{n}.filtered.faa`) from each MAG as input and performs two `DIAMOND` searches:

1. **Search against the NR database** to identify homologous sequences. The results are saved as:

   ```text
   max300_subseq.fmt6
   ```

2. **Self-search against the corresponding MAG protein database** to identify within-genome homologs. The results are saved as:

   ```text
   {n}.fmt6
   ```

where `{n}` is the name of the corresponding MAG.

### Generate Taxonomic Lineages

Extract the taxonomic IDs (the last column) from `max300_subseq.fmt6` and use `TaxonKit` to retrieve the corresponding taxonomic lineages.

```text
echo 112416 | taxonkit lineage | taxonkit reformat | cut -f 1,3
```

The resulting lineage information should be saved as `reformat.txt` in the following format:

```text
112416    Eukaryota;Ascomycota;Lecanoromycetes;Lecanorales;Parmeliaceae;Letharia;Letharia columbiana
```

where the first column is the NCBI Taxonomy ID and the second column is the corresponding taxonomic lineage from kingdom to species.

### AI calculation

Use the two alignment result files generated in the previous steps together with the lineage file as input to run `xxx.py`.

```text
python3 ai_calc.py {n}.fmt6 max300_subseq.fmt6 reformat.txt
```

This script generates an `{n}_ai.txt` file for each MAG. This file contains the Alien Index (AI) calculation results and related information for all genes in the corresponding sample.

The script also extracts genes that meet the predefined filtering criteria together with their homologous sequences, generating one `{prot}_prealign.fa` file for each candidate gene. These FASTA files are used as input for subsequent sequence alignment and phylogenetic analyses.

## Step 3. Tree Inference

Before tree construction, the candidate protein sequences require additional processing.

For each `{prot}_prealign.fa` file, the following steps are performed:

1. Align the sequences using `MAFFT`.
2. Trim the resulting alignment using `trimAl`.
3. Construct a phylogenetic tree using `IQ-TREE`.

These steps can be executed by running the provided script, for example:

```text
python3 maf_trimal_iqtree.py prealign/ --step mafft
```

For each `IQ-TREE` output (`*.treefile`), reroot the tree using the `get_midpoint_outgroup()` function provided by the Python package `ete3`, and save the rerooted tree as a new Newick file.

Upload the rerooted trees to `iTOL` for manual inspection. The phylogenetic topology is then evaluated to determine whether each candidate HGT gene should be retained for downstream analyses.

## Step 4. Confidence Classification

### Input Preparation

At this stage, a set of candidate HGT genes has been identified.

First, collect the protein sequences of all candidate HGT genes into a single FASTA file named `hgt.fa`.

Next, obtain all taxonomy IDs belonging to the taxonomic subtree of the target lineage using `TaxonKit`. For example, to analyze HGT in fungi, run:

```bash
taxonkit list --ids 4751
```

Remove any unnecessary whitespace from the output and save the resulting taxonomy IDs to a file named `taxid`. The file should contain one taxonomy ID per line, for example:

```text
4751
57731
42900
45238
84418
84419
84420
84421
84422
84423
```

### Generate the Gene–Contig Mapping File

Run `stat_hgtincontigs.py` using the input files prepared in the previous steps.

```bash
python3 stat_hgtincontigs.py
```

The script generates a file named `gene_on_contig`. The `gene_on_contig0422` file used in this study has been included in this repository.

## Step 5. Donor Inference

### Input Preparation

Before proceeding to the next step, prepare the following two input files.

Generate a file named `treefile_hgt` containing only the phylogenetic trees of the candidate HGT genes. The file format should be identical to the `treefile` generated in **Step 1**.

Generate a file named `hqmag.taxon` containing the MAG accession and the corresponding taxonomic lineage. The file format should be the same as `reformat.txt` and can be generated using `TaxonKit`.

An example entry is shown below:

```text
GCA_964255575    Eukaryota;Ascomycota;Lecanoromycetes;Peltigerales;Collemataceae;Leptogium;Leptogium austroamericanum
```

The first column contains the MAG accession, and the second column contains the corresponding taxonomic lineage.

### Generate Donor file

After preparing all required input files, run:

```bash
python3 donor_species.py
```

The script uses the following input files generated in the previous steps:

- The rerooted phylogenetic trees.
- The `DIAMOND` NR search results (`max300_subseq.fmt6`).
- The `DIAMOND` self-search results (`{n}.fmt6`).
- The `{n}_ai.txt` files.
- The taxonomic lineage file (`reformat.txt`).
- Other intermediate files generated during the pipeline.

> **Note:** Ensure that all input file paths are correctly configured before running the script. Also, the current implementation of `donor_species.py` is configured for fungal HGT analyses. For photobiont HGT analyses, modify the `find_top_candidate` function by replacing `Viridiplantae` with `Fungi` so that the target lineage is correctly identified.

Among the output files, `donor_species_final` is the final result used for all downstream analyses. The `donor_species_total.xlsx` file included in this repository is the donor annotation file used  in this study. 

## Step 6. Bioinformatic validation of HGT-acquired genes

Place the raw sequencing FASTQ files in the `sample_fq` directory. The other required input file, `hq.name`, was generated in a previous step.

After preparing the input files, run:

```bash
python3 breakpoint_bwa_samtools.py
```

> **Note:** The `cmd` string in the script is configured for our laboratory's computing environment and may not run successfully on other systems. If you wish to reproduce the analysis, please modify the command according to your local environment while keeping the same software, parameters, and workflow.




