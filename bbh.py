#!/usr/bin/env python3
"""
@author: Lucas De Vrieze (r0665032)
Comparative genomics - graded assignment

Orthology, paralogy and co-orthology
"""

# See PSB file blast_both_all.pbs for bash commands to generate the reciprocal blast results files.

import csv
from collections import Counter

## Reading result files
with open("TCqGVd_a.out", "r") as TCGV_handle:
    TCGV_reader = csv.reader(TCGV_handle, delimiter="\t")
    TCGV = [entry for line in TCGV_reader for entry in line]
with open("GVqTCd_a.out", "r") as GVTC_handle:
    GVTC_reader = csv.reader(GVTC_handle, delimiter="\t")
    GVTC = [entry for line in GVTC_reader for entry in line]
with open("GVqGVd_a.out", "r") as GVGV_handle:
    GVGV_reader = csv.reader(GVGV_handle, delimiter="\t")
    GVGV = [entry for line in GVGV_reader for entry in line]
with open("TCqTCd_a.out", "r") as TCTC_handle:
    TCTC_reader = csv.reader(TCTC_handle, delimiter="\t")
    TCTC = [entry for line in TCTC_reader for entry in line]

## Find orthologs via BBHs
# getting lists of best hits for each protein of both organisms from BLAST output
TC_besthits = dict(zip(reversed(TCGV[0::12]), reversed(TCGV[1::12])))
GV_besthits = dict(zip(reversed(GVTC[0::12]), reversed(GVTC[1::12])))
# combine to find BBHs
orthologs = {}
for h in GV_besthits.keys():
    try:
        if TC_besthits[GV_besthits[h]] == h:
            orthologs[h] = GV_besthits[h]
    except:
        continue
    
# nr. of protein coding genes: count number of entries in Fasta files
# grep -c '>' GloeobacterViolaceus.faa -> 4452
# grep -c '>' TheobromaCacao.faa -> 30854
# Some proteins didn't generate sufficiently strong hits (GV: 4392 best hits; TC: 30599 best hits)
# nr. of orthologs = nr. of BBHs -> 1991
        
## Find paralogs via significant self-BLASTing hits
# Getting lists of hit self-pairs with probabilities
TCTC_probs = dict(zip(list(zip(TCTC[0::12], TCTC[1::12])), TCTC[10::12]))
GVGV_probs = dict(zip(list(zip(GVGV[0::12], GVGV[1::12])), GVGV[10::12]))
# Extracting significant self-pairs to get paralogs
TC_paralogs = [pair for pair,prob in TCTC_probs.items() if float(prob) <= 0.00001 and pair[0] != pair[1]]
TC_prots_with_paralogs = list(set([pair[0] for pair in TC_paralogs]))
GV_paralogs = [pair for pair,prob in GVGV_probs.items() if float(prob) <= 0.00001 and pair[0] != pair[1]]
GV_prots_with_paralogs = list(set([pair[0] for pair in GV_paralogs]))
        
## Proteins in TC homologous with a GV protein
# Getting all TC-GV pairs with probabilities from BLAST output
TCGV_probs = dict(zip(list(zip(TCGV[0::12], TCGV[1::12])), TCGV[10::12]))
GVTC_probs = dict(zip(list(zip(GVTC[0::12], GVTC[1::12])), GVTC[10::12]))

# Selecting groups of TC proteins significantly co-orthologous with a GV protein, that have a BBH and at least one paralog
# Filtering for groups not larger than 5
counts_all_coorthologs = dict(Counter(TC_besthits.values()))
coorthologs = {}
for k in counts_all_coorthologs.keys():
    if k in orthologs.keys():
        group = [pair for pair,prob in TCGV_probs.items() if k == pair[1] and pair[0] in TC_prots_with_paralogs and float(prob) <= 0.00001]
        if len(group) > 1 and len(group) <= 5:
            coorthologs[k] = group

# Just picking one from the list: WP_011140090.1
# Find the TC proteins that point to WP_011140090.1 as their best hit
GV_prot = "WP_011140090.1"
TC_prot = [p[0] for p in coorthologs[GV_prot]]
prot_ortholog = orthologs["WP_011140090.1"]
print(GV_prot)
print("Co-orthologs:\t" + str(TC_prot))
print("BBH:\t" + prot_ortholog)

# Writing orthologs results file
with open("orthologs.txt", "w") as res_file:
    for GV,TC in orthologs.items():
        res_file.write(GV + "\t" + TC + "\n")
            
# Writing paralogs results file
with open("paralogs.txt", "w") as res_file:
    for pair in TC_paralogs:
        res_file.write(pair[0] + "\t" + pair[1] + "\n")
    for pair in GV_paralogs:
        res_file.write(pair[0] + "\t" + pair[1] + "\n")
        
# Writing co-orthologs results file
with open("coorthologs.txt", "w") as res_file:
    for pairs in coorthologs.values():
        for pair in pairs:
            res_file.write(pair[0] + "\t" + pair[1] + "\n")