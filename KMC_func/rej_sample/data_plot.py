# plot generated data from rejection sampling cpp using Boltzmann factor PDF

import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import glob, os
import seaborn as sns

parentroot = "./"
os.chdir(parentroot)

colnames = ["s", "pdf"]
colnames2 = ["data", "norm",""]

filelist = glob.glob("*.txt")
filelist = list(set([x[:-8] for x in filelist]))

for prefix in filelist:
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
    hist, bins = np.histogram(samples, bins='auto',density=True)
    plt.bar(bins[:-1], hist * normalFactor, width=np.diff(bins), align='edge')


    fig.savefig(f"{prefix}.jpeg", dpi=100)

    print(prefix)