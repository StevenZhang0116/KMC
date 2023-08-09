# run this file in KMC_func directory
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import glob, os
from mpl_toolkits import mplot3d
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
import matplotlib.tri as tri
import math

parentroot = "../build/"
os.chdir(parentroot)

file = "savedata_order6.txt"
# file2 = "savedata.txt"
colnames = ["id","alpha","freelength","err","time","space"]
colnames2 = ["alpha","freelength","err","time"]

df = pd.read_csv(file, sep=",", on_bad_lines='skip', names=colnames2, header=None)
# df2 = pd.read_csv(file2, sep=",", on_bad_lines='skip', names=colnames2, header=None)

# figure 1
fig, ax = plt.subplots()
plt.scatter(df["alpha"], df["freelength"], c = df["err"], cmap = "jet", norm=colors.LogNorm())
cbar = plt.colorbar(label='error')
# cbar.ax.set_yscale('log')
# cbar.ax.set_ylabel('Absolute Error Ratio (log scale)')
# cbar.mappable.set_clim(vmin=pow(10,-15), vmax=pow(10,1))
ax.set_xlabel("alpha (mu m^(-2))")
ax.set_ylabel("freelength (mu m)")

fig.savefig("data-order-6-error.jpeg", dpi=100)
plt.clf()

# figure 2
fig, ax = plt.subplots()
plt.scatter(df["alpha"], df["freelength"], c = df["time"], cmap = "jet", norm=colors.LogNorm())
cbar = plt.colorbar(label='error')
cbar.ax.set_ylabel('Build Time Ratio (log scale)')
# cbar.mappable.set_clim(vmin=pow(10,-4), vmax=pow(10,1))
ax.set_xlabel("alpha (mu m^(-2))")
ax.set_ylabel("freelength (mu m)")

fig.savefig("data-order-6-time.jpeg", dpi=100)
plt.clf()