import sys
import os
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd

matplotlib.rc('font', size=13)
#plt.rcParams['figure.constrained_layout.use'] = True
sns.set_style('white')

files = sys.argv[1:]

def parse_time(f):

    d = {}
    #[echtvar] evaluated 3533 variants (1542 / second). wrote 6 variants.

    for line in open(f):
        if 'Elapsed' in line:
            d['wall time orig'] = line.split('):')[1].strip().split(':')
            wt = 0.0
            for i, v in enumerate(reversed(d['wall time orig'])):
               wt += max(1, 60 * i) * float(v)
            d['milliseconds'] = 1000 * wt
            d.pop('wall time orig')

        """
        if 'Maximum resident' in line:
            d['Memory (MB)'] = int(line.split('):')[1].strip()) / 1000.0
        if 'User time' in line:
            d['User time (seconds)'] = float(line.split('):')[1].strip())
        if line.startswith('#SIZE:'):
            d['File size (GB)'] = float(line.split()[1]) / 1000.0
        """
        if line.startswith("[echtvar]"):
            d['variants'] = int(line.split(". wrote ")[1].split()[0])
    return d

import numpy as np

# https://stackoverflow.com/a/51535326
def show_values_on_bars(axs):
    def _show_on_single_plot(ax):        
        for p in ax.patches:
            _x = p.get_x() + p.get_width() / 2
            _y = p.get_y() + p.get_height()
            value = '{:.1f}'.format(p.get_height())
            ax.text(_x, _y, value, ha="center") 

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _show_on_single_plot(ax)
    else:
        _show_on_single_plot(axs)


tools = [parse_time(f) for f in files]



df = pd.DataFrame.from_records(tools)
df.to_csv('icgc.plots.tsv', index=False, sep="\t")



fig, axes = plt.subplots(1, 2, figsize=(8, 4)) #, constrained_layout=True)

axes[0].hist(df.milliseconds / 1000.0, 50)
axes[0].set_xlabel("Echtvar run-time (seconds)")
axes[0].set_ylabel("Count")

print(df.shape)
print(df[df.variants > 100])
df = df[df.variants < 100]
axes[1].hist(df.variants, 50)
axes[1].set_xlabel("Number of variants after filtering")
axes[1].set_ylabel("Count")

plt.tight_layout()
plt.savefig('echtvar-somatic-filter.png', format='png', dpi=1200)
plt.savefig('echtvar-somatic-filter.eps', format='eps', dpi=1200)
plt.show()
