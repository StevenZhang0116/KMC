import glob, os
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.pyplot import figure

figure(figsize=(8, 6), dpi=80)

parentroot = "../build/int-res/"
outputname = "test.png"
try: 
    os.remove(parentroot+outputname)
except: 
    pass

def calAlpha(x, i):
    res = x / (2 * 0.00411)
    cc = len(str(int(res)))
    return round(x / (2 * 0.00411), 6-cc)

os.chdir(parentroot)
files = []
for file in glob.glob("*.txt"):
    if file != "CMakeCache.txt": 
        # trasf = float(format(float(file[:-4]),'.6f'))
        files.append(file[:-4])
# files.sort()


colnames = ['integral','distance','ver-dist','alpha']
allpd = pd.DataFrame(columns=colnames)

for i in range(len(files)):
    file = files[i]
    # file = format(float(file),'.6f') + ".txt"
    file = file+".txt"
    leg = file[:-4]
    df = pd.read_csv(file, sep=",", on_bad_lines='skip', names=colnames, header=None)
    frames = [allpd, df]
    allpd = pd.concat(frames)

unique_vertical = allpd['ver-dist'].unique()

for ival in unique_vertical:
    subdf = allpd[allpd['ver-dist'] == ival]
    plt.plot(subdf['distance'],subdf['integral'],marker="o", linestyle="-",linewidth=3)


plt.savefig(outputname)
