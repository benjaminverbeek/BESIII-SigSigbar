# Simple programme to read root file

# Using uproot
import uproot
import numpy as np
import matplotlib.pyplot as plt

# Save figures here
saveLoc = r"C:\Users\benja\OneDrive\Documents\15c project\USTC\Ana_plots_rhopi"

print(f"{f' Plotting histograms for rhopi_ana.root ':=^60}")
# Open root file
file = uproot.open(r"C:\Users\benja\OneDrive\Documents\15c project\USTC\RhopiAlg\TestRelease\TestRelease-00-00-95\run\rhopi_ana.root")

count = 0
print('...')

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
        print(f"Plotting data: {branchName} in {treeName}... \t\t               ", end='\r')
        # Plot data as histogram without filled bins
        plt.hist(data, bins=100, histtype='step')
        plt.title("Histogram of branch: " + branchName + " in tree: " + treeName)
        plt.xlabel("Data")
        plt.ylabel("Count")

        # save figure as png to C:\Users\benja\OneDrive\Documents\15c project\USTC\Ana_plots_rhopi
        plt.savefig(saveLoc + r'\hist_' + treeName + '_' + branchName + ".png")
        count += 1
        # clear hist for next iteration
        plt.clf()

        # Print and overwrite same line
        print(f"Plotting data: {branchName} in {treeName}... \t\t Done!", end='\r')

# center text surrounded by '='
print(f"{f' Plotted {count} histograms! ':=^60}")

