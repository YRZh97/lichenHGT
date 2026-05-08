# lichenHGT
**codes for HGT detection in lichen (Crude code, incomplete, use with caution)**

kmerfreq_exclud.py ------- calculate k-mer frequency of each contig and exclude outliers, corresponding to 'Contig filtering and contamination detection' part in Fig.2  
process_fcs_res.py ------- flag contaminant sequences, corresponding to 'Contig filtering and contamination detection' part in Fig.2  
filt.py ------------------ select sequences based on the results of the previous two scripts, corresponding to 'Contig filtering and contamination detection' part in Fig.2  
diamond.py                   sequence alignment, corresponding to 'Alien Index score calculation' part in Fig.2  
ai_calc.py                  AI calculation, corresponding to 'Alien Index score calculation' part in Fig.2  
maf_trimal_iqtree.py        tree construction based on AI screening, corresponding to 'Tree inference' part in Fig.2  
stat_hgtincontigs.py        assign confidence levels to candidate genes, corresponding to 'Confidence classification' part in Fig.2  
donor_species.py            donor inference, corresponding to 'Donor inference' part in Fig.2  
breakpoint_bwa_samtools.py  obtain coverage information at the breakpoint, corresponding to Fig. 5b &5c  
depth.py                    coverage calculation, corresponding to Fig. 5b &5c
