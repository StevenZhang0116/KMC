# plot generated data from rejection sampling cpp using Boltzmann factor PDF

import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import glob, os
import seaborn as sns

parentroot = "./data/"
os.chdir(parentroot)

colnames = ["s", "pdf"]
colnames2 = ["data", "norm","totaltime","dt"]

filelist = glob.glob("*.txt")
filelist = list(set([x[:-8] for x in filelist]))

index = 1

for prefix in filelist:
    if index == 0: 
        # handle real data
        df = pd.read_csv(f"{prefix}data.txt", sep=",", on_bad_lines='skip', names=colnames, header=None)
        fig = plt.figure()
        plt.plot(df["s"],df["pdf"], linestyle='-.', color='red')
        plt.xlabel("s")
        plt.ylabel("integral)")
        plt.title("PDF")

        # # handle resampled data
        # df2 = pd.read_csv(f"{prefix}samp.txt", sep=",", on_bad_lines='skip', names=colnames2, header=None)
        # samples = df2["data"]
        # normalFactor = df2["norm"][0]
        # ax = sns.histplot(samples, kde=True, stat="density",element="step", fill=False)
        # hist_data = ax.get_lines()[0].get_data()
        # hist_data_y_normalized = hist_data[1] * normalFactor
        # plt.plot(hist_data[0], hist_data[1], color='blue')
        # plt.title('Histogram with Normalized Frequency')

        df2 = pd.read_csv(f"{prefix}samp.txt", sep=",", on_bad_lines='skip', names=colnames2, header=None)
        samples = df2["data"]
        normalFactor = df2["norm"][0]
        hist_values, bin_edges = np.histogram(samples, bins='auto',density=True)
        normalizedHists = hist_values * normalFactor
        plt.bar(bin_edges[:-1], normalizedHists, width=np.diff(bin_edges))

        suffix = "-resampled"

    elif index == 1:
        df = pd.read_csv(f"{prefix}samp.txt", sep=",", on_bad_lines='skip', names=colnames2, header=None)
        fig = plt.figure()
        data = df["data"]
        totalTime = int(df["totaltime"][0])
        dt = int(df["dt"][0])
        timeInterval = int(totalTime / dt)
        totalSample = len(data)
        grid = np.linspace(0, totalTime, totalSample)
        samplePerGrid = int(totalSample / timeInterval)

        bind = 0
        unbind = 0
        for i in range(timeInterval):
            aslice = data[i * samplePerGrid: (i + 1) * samplePerGrid]
            if aslice.iloc[-1] == 1:
                bind += 1
            else:
                unbind += 1

        print(f"Bind/Unbind Ratio: {bind/unbind}")
       
        plt.plot(grid, data)
        plt.title(f"Total Sample = {totalSample}; Total Time = {totalTime}; Interval = {dt}")

        suffix = "-mc"

    fig.savefig(f"{prefix}{suffix}.jpeg", dpi=100)
