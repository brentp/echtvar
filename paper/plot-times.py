import sys
import os
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd

sns.set_style('white')

files = sys.argv[1:]

def parse_time(f):

    d = {'tool': os.path.basename(f).split('-')[0]}

    for line in open(f):
        if 'Elapsed' in line:
            d['wall time orig'] = line.split('):')[1].strip().split(':')
            wt = 0.0
            for i, v in enumerate(reversed(d['wall time orig'])):
               wt += max(1, 60 * i) * float(v)
            d['Wall time (seconds)'] = wt
            d.pop('wall time orig')

        if 'Maximum resident' in line:
            d['Memory (MB)'] = int(line.split('):')[1].strip()) / 1000.0
        if 'User time' in line:
            d['User time (seconds)'] = float(line.split('):')[1].strip())
    return d


tools = [parse_time(f) for f in files]
df = pd.DataFrame.from_records(tools)

df = pd.melt(df, id_vars=('tool',),
    value_vars=('Wall time (seconds)', 'User time (seconds)', 'Memory (MB)'), value_name=('value'))

fig, axes = plt.subplots(3, 1, figsize=(4, 12))

for i, var in enumerate(('Wall time (seconds)', 'User time (seconds)', 'Memory (MB)')):

    sub = df.loc[df.variable == var, :]
    print(sub.head())
    sns.barplot(data=sub, x='tool', y='value', ax=axes[i])
    axes[i].set_ylabel(var)
    axes[i].set_xlabel(None)

plt.tight_layout()
plt.show()
