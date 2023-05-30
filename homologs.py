#!/usr/bin/env python3
"""
@author: Lucas De Vrieze (r0665032)
Comparative genomics - graded assignment

Getting sequences of homologs in 25 handpicked proteomes
"""

# See PSB file get_homologs.pbs for bash commands to detect the homologs in all 25 species by blast'ing

from Bio import SeqIO
import os

with open("homologs_ids", "r") as id_handle:
    ids = id_handle.read().split("\n")[:-1]
    
# Scan all proteome fasta files for the entries to get, make a user-friendly label and export them to a new fasta file
proteomes = os.listdir('proteomes')
records = []
for fasta in proteomes:
    proteome = SeqIO.to_dict(SeqIO.parse('/'.join(["proteomes",fasta]), "fasta"))
    for ID in ids:
        if ID in proteome.keys():
            orig = proteome[ID]
            new_id = fasta.split('.')[0] + ':' + ID
            # Keeping out a partial sequence without start codon
            if orig.seq[0] != 'M':
                continue
            record = SeqIO.SeqRecord(orig.seq, id = new_id, name = new_id, description = orig.description)
            records.append(record)
with open("homologs.fasta", "w") as handle:
    SeqIO.write(records, handle, "fasta")