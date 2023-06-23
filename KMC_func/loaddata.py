import glob, os
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 

parentroot = "../build/int-res/"
outputname = "test.png"
try: 
    os.remove(outputname)
except: 
    pass

os.chdir(parentroot)
plt.figure()
files = []
for file in glob.glob("*.txt"):
    if file != "CMakeCache.txt": 
        trasf = float(format(float(file[:-4]),'.6f'))
        files.append(trasf)
files.sort()

for file in files:
    colnames = ["integral", "sbound"]
    file = format(float(file),'.6f') + ".txt"
    leg = file[:-4]
    df = pd.read_csv(file, sep=",", on_bad_lines='skip', names=colnames, header=None)
    sbound = df["sbound"]
    integral = df["integral"]
    plt.plot(sbound, integral, linewidth=2, label=str(leg))

plt.legend()
plt.xlabel("sbound")
plt.ylabel("integral")
plt.savefig(outputname)
        
    