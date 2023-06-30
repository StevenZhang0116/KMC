import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import glob, os
from mpl_toolkits import mplot3d

parentroot = "../build/3d-data"
os.chdir(parentroot)

colnames = ["x","y","val"]

fig = plt.figure()
fig.set_size_inches(18.5, 10.5, forward=True)
ax = plt.axes(projection='3d')

for file in glob.glob("*.txt"):
    oname = file[:-4]
    df = pd.read_csv(file, sep=",", on_bad_lines='skip', names=colnames, header=None)

    ax = plt.axes(projection='3d')
    ax.scatter(df["x"], df["y"], df["val"], c='b', marker='o');

    ax.set_xticks(np.round(np.linspace(min(df["x"]), max(df["x"]), 10)))
    ax.set_yticks(np.round(np.linspace(min(df["y"]), max(df["y"]), 10)))
    ax.set_zticks(np.round(np.linspace(min(df["val"]), max(df["val"]), 10), 4))

    ax.set_xlabel("Vertical Distance", fontsize = 16)
    ax.set_ylabel("Scan Length", fontsize = 16)
    ax.set_zlabel("Lookup Table Result", fontsize = 16)
    ax.set_title(f"{oname}", fontsize = 16)

    fig.savefig(f"{oname}-3D-result.jpeg", dpi=100)

    