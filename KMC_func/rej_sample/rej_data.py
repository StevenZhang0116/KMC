# plot generated data from rejection sampling cpp examples

import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
import seaborn as sns
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

dimension_index = 1

if dimension_index == 1: 
    fig = plt.figure()
    colnames = ["data",""]
    prefix = "rej_sample_data_1d"
    df = pd.read_csv(f"{prefix}.txt", sep=",", on_bad_lines='skip', names=colnames, header=None)
    samples = df["data"]

    def target_distribution(x):
        mean = 3.0
        variance = 3.0
        constant = 0.0
        return np.exp(-0.5 * ((x - mean) / variance) * ((x - mean) / variance)) / (variance * np.sqrt(2 * np.pi)) + constant

    num = 10000
    xkk = np.linspace(-10, 15, num)
    ykk = [target_distribution(x) for x in xkk]
    # plot 1D normal distribution
    plt.plot(xkk,ykk,color="red")

    # try to plot approximated 1D normal distribution
    sns.histplot(samples, kde=True, stat="density")

    fig.savefig(f"{prefix}.jpeg", dpi=100)

elif dimension_index == 2:
    fig = plt.figure()
    colnames = ["data1","data2",""]
    prefix = "rej_sample_data_2d"

    df = pd.read_csv(f"{prefix}.txt", sep=",", on_bad_lines='skip', names=colnames, header=None)
    samples_x = df["data1"]
    samples_y = df["data2"]
    
    def target_distribution(x, y, mean_x, mean_y, var_x, var_y):
        exponent = -0.5 * ((np.power((x - mean_x) / var_x, 2)) + (np.power((y - mean_y) / var_y, 2)))
        normalization = 1.0 / (2.0 * np.pi * var_x * var_y)
        return normalization * np.exp(exponent)

    mean_x = 0.0
    mean_y = 0.0
    var_x = 2.0
    var_y = 2.0

    x_min = -10.0; 
    x_max = 10.0;  
    y_min = -10.0; 
    y_max = 10.0;  
    num = 100

    xkk = np.linspace(x_min,x_max,num)
    ykk = np.linspace(y_min,y_max,num)

    xres = [x for x in xkk for y in ykk]
    yres = [y for x in xkk for y in ykk]
    zres = [target_distribution(x, y, mean_x, mean_y, var_x, var_y) for x in xkk for y in ykk]

    # # === 3D Showcase === 
    # fig = plt.figure(figsize=(12, 5))
    # ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    # ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    # # plot 2D normal distribution
    # ax1.scatter3D(xres, yres, zres, color = "blue")

    # # try to plot approximated 2D normal distribution
    # hist, xedges, yedges = np.histogram2d(samples_x, samples_y)
    # frequencies = hist.ravel()
    # frequencies /= frequencies.sum()
    # X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
    # dx = (xedges[1] - xedges[0])
    # dy = (yedges[1] - yedges[0])
    # ax2.bar3d(X.ravel(), Y.ravel(), np.zeros_like(frequencies), dx, dy, frequencies)


    # === 2D Showcase ===
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    # plot sample data
    sns.histplot(x=samples_x,y=samples_y,kde=True,ax=axs[0],cmap="Blues",cbar=True,stat="density",bins=100)
    # sns.kdeplot(x=samples_x,y=samples_y,cbar=True,cmap="Blues",fill=True,ax=axs[0])
    axs[0].set_xlim(-10,10)
    axs[0].set_ylim(-10,10)
    axs[0].set_title("Sampled Data", fontsize=14)

    # plot target distribution data
    sns.histplot(x=xres,y=yres,weights=zres,ax=axs[1],cbar=True,cmap="Blues",stat="density",bins=100)
    # axs[1].scatter(xres, yres, c=zres, cmap="Blues", cbar=True)/
    axs[1].set_xlim(-10,10)
    axs[1].set_ylim(-10,10)
    axs[1].set_xlabel("data1", fontsize=12)
    axs[1].set_ylabel("data2", fontsize=12)
    axs[1].set_title("Targeted Distribution", fontsize=14)

    fig.savefig(f"{prefix}.jpeg", dpi=100)