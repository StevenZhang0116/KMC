# run this file in KMC_func directory
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import glob, os
from mpl_toolkits import mplot3d
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D

parentroot = "../build/3d-data"

# del current images
for filename in os.listdir(parentroot):
    if filename.endswith('.jpeg') or filename.endswith('.jpg'):
        file_path = os.path.join(parentroot, filename)
        os.remove(file_path)
        print(f"Removed file: {file_path}")


os.chdir(parentroot)

colnames = ["x","y","val","error","integral"]

fig = plt.figure()
# fig.set_size_inches(18.5, 10.5, forward=True)

# Create custom LogFormatter
class PowerLogFormatter(ticker.LogFormatter):
    def __call__(self, x, pos=None):
        return f'$10^{{{int(np.log10(x))}}}$'

for file in glob.glob("*.txt"):
    oname = file[:-4]
    df = pd.read_csv(file, sep=",", on_bad_lines='skip', names=colnames, header=None)

    # plot 3D graph
    # ax = plt.axes(projection='3d')
    # scatter = ax.scatter(df["x"], df["y"], df["val"], c=df["val"], marker='o',vmin=min(df["val"]),vmax=max(df["val"]));
    # cbar = fig.colorbar(scatter)

    # ax.set_xticks(np.round(np.linspace(min(df["x"]), max(df["x"]), 10)))
    # ax.set_yticks(np.round(np.linspace(min(df["y"]), max(df["y"]), 10)))
    # ax.set_zticks(np.round(np.linspace(min(df["val"]), max(df["val"]), 10), 2))

    # ax.set_xlabel("Vertical Distance", fontsize = 12)
    # ax.set_ylabel("Scan Length", fontsize = 12)
    # ax.set_zlabel("Lookup Table Result", fontsize = 12)
    # # ax.set_title(f"{oname}", fontsize = 12)

    # fig.savefig(f"{oname}-3D-result.jpeg", dpi=100)
    # plt.clf()

    # error plot
    # fig, ax = plt.subplots()
    # plt.scatter(df["x"], df["y"], c = df["error"], cmap = "jet", norm=colors.LogNorm())
    # cbar = plt.colorbar(label='error')
    # cbar.ax.set_yscale('log')
    # cbar.ax.set_ylabel('Error (log scale)')
    # cbar.mappable.set_clim(vmin=pow(10,-15), vmax=pow(10,1))

    # plt.axis('equal')
    # ax.set_title(f"{oname}", fontsize = 12)
    # ax.set_xlabel("r_{perp} - Perpendicular Distance", fontsize=12)
    # ax.set_ylabel("s - Scan Length",fontsize=12)
    # fig.savefig(f"{oname}-error-contour.jpeg", dpi=100)
    # plt.clf()

    # # boundary plot
    # fig, ax = plt.subplots()
    # threshold = 0.00
    # # filter everything smaller than threshold to 0
    # df["val"] = df["val"].apply(lambda x: 0 if x < threshold else x)
    # # scatter
    # scatter = ax.scatter(df["x"], df["y"], c = df["val"], cmap = 'viridis')
    # # scatter = ax.scatter(df["x"], df["y"], df["val"], c=df["val"], cmap='viridis', marker='o');
    # cbar = plt.colorbar(scatter)
    # fig.savefig(f"{oname}-boundary-result.jpeg", dpi=100)
    # plt.clf()

    # # simple plot
    # fig, ax = plt.subplots()
    # scatter = ax.scatter(df["x"], df["y"], c = df["val"], cmap = 'viridis')
    # fig.savefig(f"{oname}-res-distribution-result.jpeg", dpi=100)

    # integral plot
    fig, ax = plt.subplots()
    plt.scatter(df["x"], df["y"], c = df["integral"], cmap = "jet")
    cbar = plt.colorbar(label='error')
    cbar.ax.set_ylabel('Error')
    cbar.mappable.set_clim(vmin=0, vmax=np.max(df["integral"]))

    plt.axis('equal')
    ax.set_title(f"{oname}", fontsize = 12)
    ax.set_xlabel("r_{perp} - Perpendicular Distance", fontsize=12)
    ax.set_ylabel("s - Scan Length",fontsize=12)
    fig.savefig(f"{oname}-integral-contour.jpeg", dpi=100)
    plt.clf()

    print(oname)



    