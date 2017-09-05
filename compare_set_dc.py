# Black box the previous version of this code
# Compare both lineages of proteins to the reference structure RMSD
# with a fixed cutoff distance

# Import the libraries
import glob
import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt

cutoff = 1.0

# Import list of names
name_set = pd.read_table('protein_list.txt', delim_whitespace =True, header=None).values
line_length = len(name_set)
lineage = 0



for lineage in range(0, 2):
    for iter in range(0, line_length):
        rmsd_filenames = glob.glob('./' + name_set[iter, lineage] + '/*.dat')
        numfiles = len(rmsd_filenames)
        # Setup the calculation of the variances of each cluster's RMSDs
        cluster_variances = np.zeros(numfiles)
        holder = pd.read_table(rmsd_filenames[0], delim_whitespace=True)
        all_rmsds = holder.iloc[:,1].values
        num_frames = len(holder)
        # Setup the calculation of the mean RMSD for that cluster
        cluster_compareRMSD = np.zeros(numfiles)
        cluster_compareRMSD[0] = np.mean(all_rmsds)
        cluster_variances[0] = np.var(all_rmsds)
        for j in range(1, numfiles):
            holder = pd.read_table(rmsd_filenames[j], delim_whitespace=True)
            holder = holder.iloc[:,1]
            all_rmsds = np.vstack((all_rmsds, holder))
            cluster_compareRMSD[j] = np.mean(all_rmsds[j,:])
            cluster_variances[j] = np.var(all_rmsds[j,:])
        del holder
            

    
    
        # Compare the RMSD from each cluster's .dat file to each other
        cluster_distribution = np.zeros(numfiles)
        cluster_AvgRMSD = np.zeros(numfiles)
        cluster_AvgNotRMSD = np.zeros(numfiles)
        cluster_notCounts = np.zeros(numfiles)
        unassigned_frames = 0.0
    
        # Take the index of the lowest value, compare to cutoff, assign if lower 
        # If not lower, unassigned frames increase
        for k in range(0, num_frames):
            index = int(all_rmsds[:,k].argmin())
            rmsd = float(min(all_rmsds[:,k]))
            if rmsd <= cutoff:
                cluster_AvgRMSD[index] += rmsd
                cluster_distribution[index] += 1
            else:
                unassigned_frames += 1
            for l in range(0, numfiles):
                if l != index:
                    cluster_AvgNotRMSD[l] += all_rmsds[l, k]
                    cluster_notCounts[l] += 1
        cluster_AvgRMSD = np.divide(cluster_AvgRMSD, cluster_distribution)
        cluster_AvgNotRMSD = np.divide(cluster_AvgNotRMSD, cluster_notCounts)
        # Append the unassigned frame count to the end of the distribution array
        cluster_distribution = np.append(cluster_distribution, unassigned_frames)
        print cluster_distribution
        print cluster_AvgRMSD
        print cluster_AvgNotRMSD
        print sum(cluster_distribution)

        np.savetxt(name_set[iter, lineage] + '_' + str(cutoff) + '_cluster_distribution.txt',cluster_distribution)
        np.savetxt(name_set[iter, lineage] + '_' + str(cutoff) + '_cluster_AvgRMSD.txt',cluster_AvgRMSD)
        np.savetxt(name_set[iter, lineage] + '_' + str(cutoff) + '_cluster_NOT_RMSD.txt', cluster_AvgNotRMSD)
        np.savetxt(name_set[iter, lineage] + '_' + str(cutoff) + '_cluster_variances.txt', cluster_variances)