import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import glob, os
from mpl_toolkits import mplot3d
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D

parentroot = "../build/3d-data"
os.chdir(parentroot)

colnames = ["x","y","val","error"]

fig = plt.figure()
fig.set_size_inches(18.5, 10.5, forward=True)
ax = plt.axes(projection='3d')

# Create custom LogFormatter
class PowerLogFormatter(ticker.LogFormatter):
    def __call__(self, x, pos=None):
        return f'$10^{{{int(np.log10(x))}}}$'


for file in glob.glob("*.txt"):
    oname = file[:-4]
    df = pd.read_csv(file, sep=",", on_bad_lines='skip', names=colnames, header=None)

    ax = plt.axes(projection='3d')
    scatter = ax.scatter(df["x"], df["y"], df["val"], c=df["val"], marker='o',vmin=min(df["val"]),vmax=max(df["val"]));
    cbar = fig.colorbar(scatter)

    ax.set_xticks(np.round(np.linspace(min(df["x"]), max(df["x"]), 10)))
    ax.set_yticks(np.round(np.linspace(min(df["y"]), max(df["y"]), 10)))
    ax.set_zticks(np.round(np.linspace(min(df["val"]), max(df["val"]), 10), 4))

    ax.set_xlabel("Vertical Distance", fontsize = 16)
    ax.set_ylabel("Scan Length", fontsize = 16)
    ax.set_zlabel("Lookup Table Result", fontsize = 16)
    ax.set_title(f"{oname}", fontsize = 16)

    fig.savefig(f"{oname}-3D-result.jpeg", dpi=100)
    plt.clf()

    plt.scatter(df["x"], df["y"], c = df["error"], cmap = "jet", norm=colors.LogNorm(vmin=min(df["error"]), vmax=max(df["error"])))
    cbar = plt.colorbar(label='error')
    cbar.ax.set_yscale('log')
    cbar.ax.set_ylabel('Error (log scale)')

    plt.axis('equal')
    fig.savefig(f"{oname}-error-contour.jpeg", dpi=100)
    plt.clf()

    print(oname)



    