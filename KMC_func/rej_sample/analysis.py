# plot generated data from rejection sampling cpp using Boltzmann factor PDF

import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import glob, os
import seaborn as sns
import dask.dataframe as dd


if __name__ == "__main__":
    parentroot = "./data/"
    os.chdir(parentroot)

    colnames = ["s", "pdf"]
    colnames2 = ["alpha","fl","data", "norm","totaltime","dt","totalsamples","theory-ratio"]
    colnames3 = ["total_samples", "time_ratio"]

    filelist = glob.glob("*.txt")
    filelist = list(set([x[:-8] for x in filelist]))

    samplepersteplist = []
    acceptratio = []

    samples_num_lst = []
    ratio_lst = []

    # plotting 
    index = 2
    for prefix in filelist:
        try: 
            if index == 0: 
                substring = "st"
                kindex = prefix.find(substring)
                oldprefix = prefix
                prefix = prefix[0: kindex-1]

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

                df2 = pd.read_csv(f"{oldprefix}samp.txt", sep=",", on_bad_lines='skip', names=colnames2, header=None)
                samples = df2["data"]
                normalFactor = df2["norm"][0]
                hist_values, bin_edges = np.histogram(samples, bins='auto',density=True)
                normalizedHists = hist_values * normalFactor
                plt.bar(bin_edges[:-1], normalizedHists, width=np.diff(bin_edges))

                suffix = "-resampled"

                print(oldprefix)
                fig.savefig(f"{oldprefix}{suffix}.jpeg", dpi=100)

            elif index == 1:
                print(prefix)
                df = pd.read_csv(f"{prefix}samp.txt", sep=",", on_bad_lines='skip', names=colnames2, header=None)
                # fig = plt.figure()
                data = df["data"]
                totalTime = int(df["totaltime"][0])
                dt = df["dt"][0]

                standard = df["theory-ratio"][0]
                alpha = df["alpha"][0]
                fl = df["fl"][0]

                timeInterval = int(totalTime / dt)
                print(timeInterval)
                totalSample = len(data)
                grid = np.linspace(0, totalTime, totalSample)
                samplePerGrid = int(totalSample / timeInterval)

                # count binding state
                bindstate = 0
                unbindstate = 0
                for i in range(timeInterval):
                    aslice = data[i * samplePerGrid: (i + 1) * samplePerGrid]
                    if aslice.iloc[-1] == 1:
                        bindstate += 1
                    else:
                        unbindstate += 1

                # count binding time (uniform grid)
                bindtime = (data == 1).sum()
                unbindtime = (data == 0).sum()

                # print(f"Bind/Unbind State Ratio: {bindstate/unbindstate}")
                # print(f"Bind/Unbind Time Ratio: {bindtime/unbindtime}")
            
                # plt.plot(grid, data)
                # plt.title(f"Total Sample = {totalSample}; Total Time = {totalTime}; Interval = {dt}")

                suffix = "-mc"

                # data loading
                substring = "st"
                kindex = prefix.find(substring) + len(substring) + 1

                samplepersteplist.append(float(prefix[kindex:]))
                acceptratio.append(bindtime/unbindtime)
            
            elif index == 2:
                df = pd.read_csv(f"{prefix}samp.txt", sep=",", on_bad_lines='skip', names=colnames3, header=None)
                samples_num = df["total_samples"][0]
                ratio = df["time_ratio"][0]
                samples_num_lst.append(samples_num)
                ratio_lst.append(ratio)
        except Exception as e:
            print(e)
            continue

    if index == 1:
        fig = plt.figure()
        sortedlist = sorted(zip(samplepersteplist,acceptratio))
        sort1, sort2 = zip(*sortedlist)
        plt.scatter(sort1,sort2)
        plt.plot(sort1,sort2)
        plt.xscale("log")
        # plt.axhline(y=standard,linestyle='--')
        plt.title(f"alpha={alpha}-fl={fl}")
        plt.xlabel("Samples per dt")
        plt.ylabel("Numerical ON/OFF Time Ratio")
        fig.savefig(f"alpha={alpha}-fl={fl}.jpeg",dpi=100)

    elif index == 2:
        fig = plt.figure()
        sortedlist = sorted(zip(samples_num_lst,ratio_lst))
        sort1, sort2 = zip(*sortedlist)
        plt.scatter(sort1,sort2)
        plt.plot(sort1,sort2)
        plt.xscale("log")
        plt.xlabel("Samples per dt")
        plt.ylabel("Numerical ON/OFF Time Ratio")
        fig.savefig(f"quick-test.jpeg",dpi=100)


        
