# Author:   Benjamin Verbeek
# Date:     2022-11-25
#           Hefei, China. USTC.
# Desription:
# Simple programme to read root file and plot the invariant masses of
# pi+, pi-, rho0 and pi0.
# It can then fit a gaussian to the data and print the mean and sigma.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Import libraries
import uproot
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

# Save figures here
saveLoc = r"C:\Users\benja\OneDrive\Documents\15c project\USTC\ana_workplots\invariant_mass"

# Base file path
baseFilePath = r"C:\Users\benja\OneDrive\Documents\15c project\USTC\all_ana_roots"

# List of files look in
filename = "rhopi_ana_1000_invMass"
plotnames = "1000_invMass"

########################################################

# MAKE PLOT FOR INVARIANT MASS

print(f"{f' Plotting histograms for {filename} ':=^80}")
# Open root file
file = uproot.open(baseFilePath + "\\" + filename + ".root")
count = 0

# fit4c
treeChoice = file.keys()[7]

# Enter chosen tree
tree = file[treeChoice]
print("Entering tree: ", tree)

# Masses of pi+, pi-, rho0 and pi0
branchChoices = [1, 2, 3, 4]
# Number of bins in histogram
nBins = [80, 10, 10, 80]

for branchChoice in branchChoices:
    # Enter chosen branch
    branch = tree[tree.keys()[branchChoice]]
    print("Entering branch: ", branch)

    # Get data
    data = branch.array()
    print(data)

    # Bin the data:
    # Get the min and max values of the data
    minVal = np.min(data)
    maxVal = np.max(data)
    # Calculate the bin width
    binWidth = (maxVal - minVal) / nBins[branchChoice - 1]
    # Make the bins
    bins = np.arange(minVal, maxVal, binWidth)
    # Bin the data
    dataHist, bins = np.histogram(data, bins=bins)
    print(dataHist)

    # Plot the data
    plt.figure(count)
    plt.plot(bins[:-1], dataHist, "r-", label="Data")
    plt.xlabel("Invariant mass (GeV)")
    plt.ylabel("Counts")
    plt.title(f"{filename} {branch.name}")
    # set xlimit as min and max values of bins
    plt.xlim(minVal, maxVal)
    #plt.show()

    # Find peaks
    peaks, _ = find_peaks(dataHist, height=10)
    # ^ indices of peaks
    print(peaks)
    pk = peaks[0] # first peak

    peakLoc = [(0, 4), (1,9), (2, 7), (2, 10)]

    # Fit a gaussian to the data, using the peaks as initial guesses
    # for the mean and sigma
    # Define the gaussian function
    # Gaussian function:
    def gaus(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    print("Peak loc: ", pk-int(0.8*pk), pk, pk+int(0.8*pk))

    # Fit the data
    popt, pcov = curve_fit(gaus, bins[peakLoc[count][0]:peakLoc[count][1]], dataHist[peakLoc[count][0]:peakLoc[count][1]], p0=[1, peaks[0], 1])
    # ^ popt is the optimal values for the parameters so that the sum of the
    # squared residuals of gaus(xdata, *popt) - ydata is minimized
    print(popt)
    # plot the fitted gaussian with a 10x denser grid
    x = np.linspace(minVal, maxVal, 10*len(bins))
    plt.plot(x, gaus(x, *popt), "b--", label="Fitted gaussian")
    plt.legend()
    plt.show()

    # branchName = str(branch).split("'")[1]
    # treeName   = str(tree).split("'")[1]
    # print(f"Plotting data: {branchName} in {treeName}...")
    # # Plot data as histogram without filled bins
    # plt.figure(branchChoice)
    # #plt.hist(data, bins=100, histtype='step')
    # plt.hist(data, bins=nBins[branchChoice-1], alpha=0.5, edgecolor='black', linewidth=1, label=plotnames)#, histtype='step')
    # plt.title("Invariant mass (branch: " + branchName + " in tree: " + treeName + ") (for rhopi)")
    # plt.xlabel("Mass [GeV]")
    # plt.ylabel("Count")
    # plt.legend()
    # # save figure as png to C:\Users\benja\OneDrive\Documents\15c project\USTC\Ana_plots_rhopi
    # plt.savefig(saveLoc + r'\hist_' + treeName + '_' + branchName + ".png")
    # plt.close()

    count += 1
    # clear hist for next iteration
    #plt.clf()

    # Print and overwrite same line
    #print(f"Plotting data: {branchName} in {treeName}... \t\t Done!", end='\r')
print("Done.")
# center text surrounded by '='
print(f"{f' Plotted {count} histograms. ':=^80}")

