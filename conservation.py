#!/usr/bin/env python3
"""
@author: Lucas De Vrieze (r0665032)
Comparative genomics - graded assignment

Conservation score metrics of the alignment of homologs from 25 species
"""

# MSA was generated using clustalO at default settings.

from Bio import AlignIO
from collections import Counter
from numpy import log2
import matplotlib.pyplot as plt
from pandas import Series

# Read alignment file
with open("homologs.msa", "r") as align_handle:
    alignment_reader = AlignIO.read(align_handle, "fasta")
    alignment = [str(entry.seq) for entry in alignment_reader]

# Calculate the conservation score
## Two-prong approach inspired by the trident score of Valdar (https://doi.org/10.1002/prot.10146)
scores = []
nr_alignments = len(alignment)
for i in range(len(alignment[0])):
    a_i = [a[i] for a in alignment]
    
    # Shannon entropy
    aa_counts = dict(Counter(a_i))
    fractions = [v/nr_alignments for v in aa_counts.values()]
    scale_factor = 1/(log2(min(nr_alignments, 21)))
    shannon = -scale_factor*sum([p*log2(p) for p in fractions])
    
    # Gappiness
    try:
        gappiness = aa_counts['-']/nr_alignments
    except:
        gappiness = 0
    
    # Combined score
    score = (1-shannon)*(1-gappiness)
    scores.append(score)

# setting up a moving average for plotting
window_size = 5
scores_series = Series(scores)
windows = scores_series.rolling(window_size)
moving_average = windows.mean().tolist()

# plotting
plt.figure()
plt.plot(moving_average)
plt.axhline(y=0.3, color='red')
plt.xlabel('Amino acid index')
plt.ylabel('Conservation score (moving average with window size 5)')
plt.xlim([0,len(alignment[0])])
plt.tight_layout()
plt.show()
plt.savefig("ConservationScore.pdf")