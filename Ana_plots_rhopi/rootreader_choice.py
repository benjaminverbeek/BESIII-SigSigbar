# Simple programme to read root file

# Using uproot
import uproot
import numpy as np
import matplotlib.pyplot as plt

# Open root file
file = uproot.open(r"C:\Users\benja\OneDrive\Documents\15c project\USTC\RhopiAlg\TestRelease\TestRelease-00-00-95\run\rhopi_ana.root")

# Print keys
[print(f'{i}: {key}', end=" | ") for i, key in enumerate(file.keys())]
#print(list(enumerate(file.keys())))
treeChoice = int(input("Which tree would you like to see?: "))

# Enter chosen tree
tree = file[file.keys()[treeChoice]]
print("Entering tree: ", tree)

# Print branches
[print(f'{i}: {key}', end=" | ") for i, key in enumerate(tree.keys())]

# Enter chosen branch
branchChoice = int(input("Which branch would you like to see?: "))
branch = tree[tree.keys()[branchChoice]]
print("Entering branch: ", branch)

# Get data
data = branch.array()
print(data)

# Plot data as histogram without filled bins
print("Plotting data...")
plt.hist(data, bins=100, histtype='step')
branchName = str(branch).split("'")[1]
treeName   = str(tree).split("'")[1]
plt.title("Histogram of branch: " + branchName + " in tree: " + treeName)
plt.xlabel("Data")
plt.ylabel("Count")

# save figure as png to C:\Users\benja\OneDrive\Documents\15c project\USTC\Ana_plots_rhopi
#plt.savefig(r"C:\Users\benja\OneDrive\Documents\15c project\USTC\Ana_plots_rhopi\hist_" + treeName + '_' + branchName + ".png")

plt.show()

