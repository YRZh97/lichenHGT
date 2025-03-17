# !/user/bin/env python
# -*- coding:utf-8 -*-

import os,sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

helpinfo = 'python3 ai_calc.py diamond905 diamondnr taxonkit'
if len(sys.argv) != 4:
    sys.exit(helpinfo)

diamond905 = sys.argv[1]
diamondnr = sys.argv[2]
taxonkit = sys.argv[3]
#seqfile = sys.argv[4]

recipient = [';Lecanoromycetes;',';Eurotiomycetes;',';Lichinomycetes;']

funginame = {}
for line in open('/data2/zhangyr/zyrtest/hgt/hqfunqi.name','r'):
    funginame['eu_' + line.strip().split(',')[0].split('.')[0]] = line.strip().split(',')[2]

taxon = {}
for line in open(taxonkit,'r'):
    if len(line.strip().split('\t')) > 1:
        taxon[line.strip().split('\t')[0]] = line.strip().split('\t')[1]

info = {}
bitscore = {}
qlength = {}

for line in open(diamond905,'r'):
    if line.strip().split('\t')[0] not in bitscore:
        bitscore.setdefault(line.strip().split('\t')[0],{line.strip().split('\t')[1]:0})
    else:
        bitscore[line.strip().split('\t')[0]][line.strip().split('\t')[1]] = 0
    if float(line.strip().split('\t')[-1]) > bitscore[line.strip().split('\t')[0]][line.strip().split('\t')[1]]:
        bitscore[line.strip().split('\t')[0]][line.strip().split('\t')[1]] = float(line.strip().split('\t')[-1])
    if diamond905.split('/')[-1].split('.')[0] + '_' + line.strip().split('\t')[0] == line.strip().split('\t')[1]:
        qlength[line.strip().split('\t')[0]] = line.strip().split('\t')[3] 

for line in open(diamondnr,'r'):
    if len(line.strip().split('\t')) == 14:
        if line.strip().split('\t')[0] not in bitscore:
            bitscore.setdefault(line.strip().split('\t')[0],{line.strip().split('\t')[1] + '|' + line.strip().split('\t')[-1].split(';')[0]:0})
        else:
            bitscore[line.strip().split('\t')[0]][line.strip().split('\t')[1] + '|' + line.strip().split('\t')[-1].split(';')[0]] = 0
        if float(line.strip().split('\t')[-3]) > bitscore[line.strip().split('\t')[0]][line.strip().split('\t')[1] + '|' + line.strip().split('\t')[-1].split(';')[0]]:
            bitscore[line.strip().split('\t')[0]][line.strip().split('\t')[1]+ '|' + line.strip().split('\t')[-1].split(';')[0]] = float(line.strip().split('\t')[-3])
    else:
        if line.strip().split('\t')[0] not in bitscore:
            bitscore.setdefault(line.strip().split('\t')[0],{line.strip().split('\t')[1] + '|': 0})
        else:
            bitscore[line.strip().split('\t')[0]][line.strip().split('\t')[1] + '|'] = 0
        if float(line.strip().split('\t')[-2]) > bitscore[line.strip().split('\t')[0]][line.strip().split('\t')[1] + '|']:
            bitscore[line.strip().split('\t')[0]][line.strip().split('\t')[1] + '|'] = float(line.strip().split('\t')[-2])
print(len(bitscore))

for query in bitscore:
    info.setdefault(query,{'maxB':0, 'bbhO':[0,'0'], 'bbhG':[0,'0'], 'fungignum':[], 'outgnum':[]})
    for subject in bitscore[query]:
        if subject.startswith('eu_') or subject.startswith('pro_'):
            if subject == diamond905.split('/')[-1].split('.')[0] + '_' + query:
                info[query]['maxB'] = bitscore[query][subject]
            else:
                if subject.rsplit('_',1)[0] in funginame:
                    info[query]['fungignum'].append(subject.rsplit('_',1)[0])
                else:
                    info[query]['outgnum'].append(subject.rsplit('_',1)[0])
                    if bitscore[query][subject] > info[query]['bbhO'][0]:
                        info[query]['bbhO'][0] = bitscore[query][subject]
                        info[query]['bbhO'][1] = subject
        else:
            if subject.split('|')[1] in taxon:
                if ';Fungi;' in taxon[subject.split('|')[1]]:
                    info[query]['fungignum'].append(subject.split('|')[1])
                    lichen = False
                    for a in recipient:
                        if a in taxon[subject.split('|')[1]]:
                            lichen = True
                            # break
                    if not lichen:
                        if bitscore[query][subject] > info[query]['bbhG'][0]:
                            info[query]['bbhG'][0] = bitscore[query][subject]
                            info[query]['bbhG'][1] = subject
                else:
                    info[query]['outgnum'].append(subject.split('|')[1])
                    if bitscore[query][subject] > info[query]['bbhO'][0]:
                        info[query]['bbhO'][0] = bitscore[query][subject]
                        info[query]['bbhO'][1] = subject
print(len(info))

tot = []
outputfile = open(diamond905.split('/')[-1].split('.')[0] + '_ai.txt','w')
outputfile.write('queryid\tqlength\tmaxB\tbbhG_subject\tbbhG_score\tbbhO_subject\tbbhO_score\treci_and_group_ids\treci_and_group_num\toutgroup_ids\toutgroup_num\talien_index\toutg_pct\n')
for key in info:
    if key in qlength:
        if len(info[key]['fungignum']) + len(info[key]['outgnum']) == 0:
            outputfile.write(key + '\t' + qlength[key] + '\t' + str(info[key]['maxB']) + '\t' + info[key]['bbhG'][1] + '\t' + str(info[key]['bbhG'][0]) + '\t' + info[key]['bbhO'][1] + '\t' + str(info[key]['bbhO'][0]) + '\t' + ';'.join(info[key]['fungignum']) + '\t' + str(len(info[key]['fungignum'])) + '\t' + ';'.join(info[key]['outgnum']) + '\t' + str(len(info[key]['outgnum'])) + '\t' + str(info[key]['bbhO'][0]/info[key]['maxB'] - info[key]['bbhG'][0]/info[key]['maxB']) + '\t0\n')   
        else:
            outputfile.write(key + '\t' + qlength[key] + '\t' + str(info[key]['maxB']) + '\t' + info[key]['bbhG'][1] + '\t' + str(info[key]['bbhG'][0]) + '\t' + info[key]['bbhO'][1] + '\t' + str(info[key]['bbhO'][0]) + '\t' + ';'.join(info[key]['fungignum']) + '\t' + str(len(info[key]['fungignum'])) + '\t' + ';'.join(info[key]['outgnum']) + '\t' + str(len(info[key]['outgnum'])) + '\t' + str(info[key]['bbhO'][0]/info[key]['maxB'] - info[key]['bbhG'][0]/info[key]['maxB']) + '\t' + str(len(info[key]['outgnum'])/(len(info[key]['fungignum']) + len(info[key]['outgnum']))) + '\n')
        if info[key]['bbhO'][0]/info[key]['maxB'] - info[key]['bbhG'][0]/info[key]['maxB'] > 0 and len(info[key]['outgnum'])/(len(info[key]['fungignum']) + len(info[key]['outgnum'])) >= 0.8:
            tot.append(key)
    else:
        pass
outputfile.close()
print(len(tot),tot)

for t in tot:
    print(t)
    lis905 = []
    for line in open(diamond905,'r'):
        if line.strip().split('\t')[0] == t:
            lis905.append(line.strip().split('\t')[1])

    gtoutput = open(t + '_prealign.fa','w')

    for record in SeqIO.parse('/data2/zhangyr/zyrtest/hgt/90_5.fas','fasta'):
        if record.id in lis905:
            SeqIO.write(record, gtoutput, 'fasta')

    for line in open(diamondnr,'r'):
        if line.strip().split('\t')[0] == t:
            if len(line.strip().split('\t')) == 14:
                taxid = line.strip().split('\t')[-1].split(';')[0]
                if taxid in taxon:
                    tax = '-'.join(taxon[taxid].split(';')[1:4]) + '-' + taxon[taxid].split(';')[-1]
                    newrecord = SeqRecord(Seq(line.strip().split('\t')[-2]), id= line.strip().split('\t')[1] + '|' + tax, description='')
                    SeqIO.write(newrecord, gtoutput, 'fasta')
    gtoutput.close()
