import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 

fig = plt.figure()
colnames = ["data",""]

df = pd.read_csv("rej_sample_data.txt", sep=",", on_bad_lines='skip', names=colnames, header=None)
plt.plot(df["data"])
fig.savefig(f"rej_sample_result.jpeg", dpi=100)