# plot generated data from rejection sampling cpp using Boltzmann factor PDF

import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import glob, os
import seaborn as sns

colnames = ["s", "pdf"]
colnames2 = ["data"]

for file in glob.glob("*.txt"):
    prefix = file[:-8]
    df = pd.read_csv(f"{prefix}data.txt", sep=",", on_bad_lines='skip', names=colnames, header=None)
    fig = plt.figure()
    plt.plot(df["s"],df["pdf"], linestyle='-.', color='red')
    plt.xlabel("s")
    plt.ylabel("integral)")
    plt.title("PDF")

    df2 = pd.read_csv(f"{prefix}samp.txt", sep=",", on_bad_lines='skip', names=colnames2, header=None)
    samples = df2["data"]
    sns.histplot(samples, kde=True, stat="density")

    fig.savefig(f"{prefix}.jpeg", dpi=100)

    print(prefix)