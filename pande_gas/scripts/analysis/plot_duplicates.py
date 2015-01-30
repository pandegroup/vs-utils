#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pp
from matplotlib import artist
import numpy as np
import sys
import seaborn as sns
sns.set(style='whitegrid')

dup = []
unq = []
with open(sys.argv[1]) as f:
  for line in f:
    dup.append(line.strip())
with open(sys.argv[2]) as f:
  for line in f:
    unq.append(line.strip())
dup = np.asarray(dup, dtype=float)
unq = np.asarray(unq, dtype=float)
print dup.size, unq.size
fig = pp.figure(figsize=(2, 2))
ax = fig.add_subplot(111)
bp = ax.boxplot([dup, unq], notch=1, positions=[0.2, 0.4], widths=0.1)
#artist.setp(bp['boxes'], facecolor=colors[i])
#artist.setp(bp['caps'], color='k')
artist.setp(bp['whiskers'], color='k', linestyle='solid')
#artist.setp(bp['medians'], color='k')
#artist.setp(bp['fliers'], color='k', marker='o', markersize=2)
ax.set_xticklabels(['Duplicate', 'Unique'])
ax.set_ylabel(r'$\Delta$ Log-Odds-Mean-AUC')
ax.set_xlim(0.1, 0.5)
ax.set_ylim(-0.5, None)
fig.savefig('duplicates.png', dpi=300, bbox_inches='tight', transparent=True)
