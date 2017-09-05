# Plot the cluster discovery and calculate the shannon entropy over the course
# of the simulations


# Importing the libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys 
import math

# Import the cluster v time dataset
cluster_v_time_filename = sys.argv[1] 
# cluster_v_time_filename = 'rnaseh_ancC_clusters_v_time.dat'
cluster_v_time = pd.read_table(cluster_v_time_filename, delim_whitespace=True)
clusterVtime = cluster_v_time.iloc[:,1]
numframes = len(clusterVtime)

# Define organism & protein name for later output files
name_split = cluster_v_time_filename.split("_")
org_prot_name = name_split[0] + '_' + name_split[1]

# Initialize relevant variables to cluster discovery
unique_clusters = np.empty(0)
cluster_discovery = np.empty(0)
unique_count = 0

# Determine uniqueness of each frame's cluster
for i in range(0, numframes):
    if clusterVtime[i] not in unique_clusters:
        unique_count += 1
        unique_clusters = np.append(arr = unique_clusters, values = clusterVtime[i])
        cluster_discovery = np.append(arr = cluster_discovery, values = unique_count)
    else: 
        cluster_discovery = np.append(arr = cluster_discovery, values = unique_count)

# Plot cluster discovery over the simulation trajectory
cluster_discovery_output_filename = org_prot_name + '_cluster_discovery.eps'
plt.plot(cluster_discovery, color = 'blue', 
         label = org_prot_name)
plt.title('Cluster discovery over the course of the simulation')
plt.xlabel('Trajectory frame (10ps resolution)')
plt.ylabel('Number of clusters')
plt.savefig(cluster_discovery_output_filename, bbox_inches='tight' )
plt.show()      

# Print the relative change in clusters found 
pct_change = float(cluster_discovery[100000])/float(cluster_discovery[124999])
print pct_change

# Remove variables that could be taking up memory space
del cluster_v_time
del clusterVtime
del pct_change

# Start of the shannon entropy calculation section of this code 
# Import the distribution over time dataset
cluster_distribution_filename = sys.argv[2]
#cluster_distribution_filename = 'rnaseh_ancC_cluster_pop.dat'
cluster_distribution_file = pd.read_table(cluster_distribution_filename, 
                                     delim_whitespace=True)
cluster_distribution = np.array(cluster_distribution_file)

# Clear file to save memory
del cluster_distribution_file

# Initialize relevant variables/arrays for shannon entropy
shannonvtime = np.empty(0)

# Calculate shannon entropy per frame, append to time evolution array
for j in range(0, numframes):
    shannon_frame = 0
    for k in range(1, unique_count + 1):
        shannon_frame += (-1)*cluster_distribution[j][k]*math.log(cluster_distribution[j][k])
    shannonvtime = np.append(arr = shannonvtime, values = shannon_frame)
    
# Plot the cluster distribution entropy
cluster_distribution_output_filename = org_prot_name + '_cluster_distribution.eps'
plt.plot(shannonvtime, color = 'blue', 
         label = org_prot_name)
plt.title('Cluster entropy over the course of the simulation')
plt.xlabel('Trajectory frame (10ps resolution)')
plt.ylabel('Entropy')
plt.savefig(cluster_distribution_output_filename, bbox_inches='tight' )
plt.show() 

pct_change = float(shannonvtime[100000])/float(shannonvtime[124999])
print pct_change
























