import glob, os
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.pyplot import figure

figure(figsize=(8, 6), dpi=80)

parentroot = "../build/rvl-res/time/"
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

# for file in files:
#     colnames = ["integral", "sbound"]
#     file = format(float(file),'.6f') + ".txt"
#     leg = file[:-4]
#     df = pd.read_csv(file, sep=",", on_bad_lines='skip', names=colnames, header=None)
#     sbound = df["sbound"]
#     integral = df["integral"]
#     plt.plot(sbound, integral, linewidth=2, label=str(leg))

# plt.legend()
# plt.xlabel("sbound")
# plt.ylabel("integral")
# plt.savefig(outputname)
        
# fig, axs = plt.subplots(2)
colnames = ['alpha','tol','time','bberr','rlerr']
allpd = pd.DataFrame(columns=colnames)

for i in range(len(files)):
    file = files[i]
    # file = format(float(file),'.6f') + ".txt"
    file = file+".txt"
    leg = file[:-4]
    df = pd.read_csv(file, sep=",", on_bad_lines='skip', names=colnames, header=None)
    allpd.loc[i] = df.iloc[0]

allpd.astype(np.uint8)
xlst = [0.1,0.5,1,5,10]
alphalist = []
for x in range(len(xlst)):
    alphalist.append(calAlpha(xlst[x],x))

for val in alphalist:
    subdf = allpd[allpd['alpha'] == val]
    subdf = subdf.sort_values(by=['tol'])
    print(subdf)
    plt.plot(subdf['tol'],subdf['time'],label=f'alpha={val}',marker="o", linestyle="-",linewidth=3)


# axs[0].plot(timepd['alpha'],timepd['time'],marker="o", linestyle="-")
# axs[1].plot(bberrpd['alpha'],bberrpd['bberr'],label='Baobzi Error',marker="o", linestyle="-")
# axs[1].plot(rlerrpd['alpha'],rlerrpd['rlerr'],label='Reverse LookUp ERROR',marker="o", linestyle="-")

# axs[1].set_xlabel('Alpha')
# axs[0].set_ylabel('Elapsed Time')
# axs[1].set_ylabel('Average Error')
# axs[1].set_yscale('log')
# axs[1].legend(loc="upper right")

plt.xlabel('Tolerance')
plt.xscale('log')
plt.ylabel('Time(s)')
plt.legend()
plt.savefig(outputname)

