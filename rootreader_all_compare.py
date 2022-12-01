# Author:   Benjamin Verbeek
# Date:     2022-11-23
#           Hefei, China. USTC.
# Desription:
# Simple programme to read root file and plot all branches.
# Can take in any number of files and plots them in the same figures.
# The purpose is to be able to modify Rhopi.cxx and see how the output
# changes.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Import libraries
import uproot
import matplotlib.pyplot as plt
import sys

# Save figures here
saveLoc = r"C:\Users\benja\OneDrive\Documents\15c project\USTC\ana_workplots"

# Base file path
baseFilePath = r"C:\Users\benja\OneDrive\Documents\15c project\USTC\all_ana_roots"

# List of files to plot
#filenames_suffix = ["1000_invMass", "100_noIPcheck_noChargeCheck", "100_noIPcheck", "100_raw", "50_raw"]
#filenames_suffix = ["1000_invMass"]
#filenames = ["rhopi_ana_" + name for name in filenames_suffix]
filenames_suffix = ["1000_4charged"]
filenames = ["sigmasigmabar_ana_" + name for name in filenames_suffix]
plotnames = filenames_suffix

# Number of bins in histogram
nBins = 50

########################################################

print(f"!!!  Programme will plot data from {len(filenames)} files  !!!")

# MAKE PLOTS
totCount = 0
for fileNr, filename in enumerate(filenames):
    print(f"{f' Plotting histograms for {filename} ':=^80}")
    # Open root file
    file = uproot.open(baseFilePath + "\\" + filename + ".root")
    count = 0

    for treeChoice in range(len(file.keys())):
        # Enter chosen tree
        tree = file[file.keys()[treeChoice]]
        #print("Entering tree: ", tree)

        for branchChoice in range(len(tree.keys())):
            # Enter chosen branch
            branch = tree[tree.keys()[branchChoice]]
            #print("Entering branch: ", branch)

            # Get data
            data = branch.array()
            #print(data)

            branchName = str(branch).split("'")[1]
            treeName   = str(tree).split("'")[1]
            sys.stdout.write("\033[K") # Clear to the end of line
            print(f"Plotting data: {branchName} in {treeName}...", end='\r')
            # Plot data as histogram without filled bins
            plt.figure(count)
            #plt.hist(data, bins=100, histtype='step')
            plt.hist(data, bins=nBins, alpha=0.5, edgecolor='black', linewidth=1)#, histtype='step')
            if fileNr == len(filenames)-1:
                plt.title("Histogram of branch: " + branchName + " in tree: " + treeName + " (for sigsigbar)")
                plt.xlabel("Data")
                plt.ylabel("Count")
                plt.legend(plotnames)
                # save figure as png to C:\Users\benja\OneDrive\Documents\15c project\USTC\Ana_plots_rhopi
                plt.savefig(saveLoc + r'\hist_' + treeName + '_' + branchName + ".png")
                plt.close()

            count += 1
            # clear hist for next iteration
            #plt.clf()

            # Print and overwrite same line
            #print(f"Plotting data: {branchName} in {treeName}... \t\t Done!", end='\r')
    sys.stdout.write("\033[K") # Clear to the end of line
    print("Done.")
    totCount += count
# center text surrounded by '='
print(f"{f' Plotted {count} histograms from a total of {totCount} datasets! ':=^80}")

