# lichenHGT
codes for HGT detection in lichen (Crude code, incomplete, use with caution)


[newworkflow.tif](https://github.com/user-attachments/files/27433871/newworkflow.tif)

workflow as follows:
kmerfreq_exclud.py ------- calculate k-mer frequency of each contig and exclude outliers, corresponding to 'Contig filtering and contamination detection' part in Fig.2
process_fcs_res.py ------- flag contaminant sequences, corresponding to 'Contig filtering and contamination detection' part in Fig.2
filt.py ------------------ sequences meeting the criteria were retained for subsequent alignment based on the results of the previous two scripts, corresponding to 'Contig filtering and contamination detection' part in Fig.2
diamond.py ---------------  sequence alignment, corresponding to 'Alien Index score calculation' part in Fig.2
ai_calc.py --------------- AI calculation, corresponding to 'Alien Index score calculation' part in Fig.2
maf_trimal_iqtree.py ----- tree construction based on AI screening, corresponding to 'Tree inference' part in Fig.2
extract_first_hit_eachgene.py ----- 
stat_hgtincontigs.py
