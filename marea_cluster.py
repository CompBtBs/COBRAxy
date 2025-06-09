# -*- coding: utf-8 -*-
"""
Created on Mon Jun 3 19:51:00 2019
@author: Narger
"""

import sys
import argparse
import os
import numpy as np
import pandas as pd
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.metrics import silhouette_samples, silhouette_score, cluster
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc   
import matplotlib.cm as cm
from typing import Optional, Dict, List

################################# process args ###############################
def process_args(args_in :List[str] = None) -> argparse.Namespace:
    """
    Processes command-line arguments.

    Args:
        args (list): List of command-line arguments.

    Returns:
        Namespace: An object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(usage = '%(prog)s [options]',
                                     description = 'process some value\'s' +
                                     ' genes to create class.')

    parser.add_argument('-ol', '--out_log', 
                        help = "Output log")
    
    parser.add_argument('-in', '--input',
                        type = str,
                        help = 'input dataset')
    
    parser.add_argument('-cy', '--cluster_type',
                        type = str,
                        choices = ['kmeans', 'dbscan', 'hierarchy'],
                        default = 'kmeans',
                        help = 'choose clustering algorythm')
    
    parser.add_argument('-k1', '--k_min', 
                        type = int,
                        default = 2,
                        help = 'choose minimun cluster number to be generated')
    
    parser.add_argument('-k2', '--k_max', 
                        type = int,
                        default = 7,
                        help = 'choose maximum cluster number to be generated')
    
    parser.add_argument('-el', '--elbow', 
                        type = str,
                        default = 'false',
                        choices = ['true', 'false'],
                        help = 'choose if you want to generate an elbow plot for kmeans')
    
    parser.add_argument('-si', '--silhouette', 
                        type = str,
                        default = 'false',
                        choices = ['true', 'false'],
                        help = 'choose if you want silhouette plots')
    
    parser.add_argument('-td', '--tool_dir',
                        type = str,
                        required = True,
                        help = 'your tool directory')
                        
    parser.add_argument('-ms', '--min_samples',
                        type = int,
                        help = 'min samples for dbscan (optional)')
                        
    parser.add_argument('-ep', '--eps',
                        type = float,
                        help = 'eps for dbscan (optional)')
                        
    parser.add_argument('-bc', '--best_cluster',
                        type = str,
                        help = 'output of best cluster tsv')
    				
    parser.add_argument(
        '-idop', '--output_path', 
        type = str,
        default='result',
        help = 'output path for maps')
    
    args_in = parser.parse_args(args_in)
    return args_in

########################### warning ###########################################
def warning(s :str) -> None:
    """
    Log a warning message to an output log file and print it to the console.

    Args:
        s (str): The warning message to be logged and printed.
    
    Returns:
      None
    """

    with open(args.out_log, 'a') as log:
        log.write(s + "\n\n")
    print(s)

########################## read dataset ######################################
def read_dataset(dataset :str) -> pd.DataFrame:
    """
    Read dataset from a CSV file and return it as a Pandas DataFrame.

    Args:
        dataset (str): the path to the dataset to convert into a DataFrame

    Returns:
        pandas.DataFrame: The dataset loaded as a Pandas DataFrame.

    Raises:
        pandas.errors.EmptyDataError: If the dataset file is empty.
        sys.exit: If the dataset file has the wrong format (e.g., fewer than 2 columns)
    """
    try:
        dataset = pd.read_csv(dataset, sep = '\t', header = 0)
    except pd.errors.EmptyDataError:
        sys.exit('Execution aborted: wrong format of dataset\n')
    if len(dataset.columns) < 2:
        sys.exit('Execution aborted: wrong format of dataset\n')
    return dataset

############################ rewrite_input ###################################
def rewrite_input(dataset :Dict) -> Dict[str, List[Optional[float]]]:
    """
    Rewrite the dataset as a dictionary of lists instead of as a dictionary of dictionaries.

    Args:
        dataset (pandas.DataFrame): The dataset to be rewritten.

    Returns:
        dict: The rewritten dataset as a dictionary of lists.
    """
    #Riscrivo il dataset come dizionario di liste, 
    #non come dizionario di dizionari
    #dataset.pop('Reactions', None)
    
    for key, val in dataset.items():
        l = []
        for i in val:
            if i == 'None':
                l.append(None)
            else:
                l.append(float(i))
   
        dataset[key] = l
    
    return dataset

############################## write to csv ##################################
def write_to_csv (dataset :pd.DataFrame, labels :List[str], name :str) -> None:
    """
    Write dataset and predicted labels to a CSV file.

    Args:
        dataset (pandas.DataFrame): The dataset to be written.
        labels (list): The predicted labels for each data point.
        name (str): The name of the output CSV file.

    Returns:
        None
    """
    #labels = predict
    predict = [x+1 for x in labels]
  
    classe = (pd.DataFrame(list(zip(dataset.index, predict)))).astype(str)

    dest = name
    classe.to_csv(dest, sep = '\t', index = False,
                      header = ['Patient_ID', 'Class'])
   
########################### trova il massimo in lista ########################
def max_index (lista :List[int]) -> int:
    """
    Find the index of the maximum value in a list.

    Args:
        lista (list): The list in which we search for the index of the maximum value.

    Returns:
        int: The index of the maximum value in the list.
    """
    best = -1
    best_index = 0
    for i in range(len(lista)):
        if lista[i] > best:
            best = lista [i]
            best_index = i
            
    return best_index
    
################################ kmeans #####################################
def kmeans (k_min: int, k_max: int, dataset: pd.DataFrame, elbow: str, silhouette: str, best_cluster: str) -> None:
    """
    Perform k-means clustering on the given dataset, which is an algorithm used to partition a dataset into groups (clusters) based on their characteristics.
    The goal is to divide the data into homogeneous groups, where the elements within each group are similar to each other and different from the elements in other groups.

    Args:
        k_min (int): The minimum number of clusters to consider.
        k_max (int): The maximum number of clusters to consider.
        dataset (pandas.DataFrame): The dataset to perform clustering on.
        elbow (str): Whether to generate an elbow plot for kmeans ('True' or 'False').
        silhouette (str): Whether to generate silhouette plots ('True' or 'False').
        best_cluster (str): The file path to save the output of the best cluster.

    Returns:
        None
    """
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    
        
    if elbow == 'true':
        elbow = True
    else:
        elbow = False
        
    if silhouette == 'true':
        silhouette = True
    else:
        silhouette = False
        
    range_n_clusters = [i for i in range(k_min, k_max+1)]
    distortions = []
    scores = []
    all_labels = []
    
    clusterer = KMeans(n_clusters=1, random_state=10)
    distortions.append(clusterer.fit(dataset).inertia_)
    
    
    for n_clusters in range_n_clusters:
        clusterer = KMeans(n_clusters=n_clusters, random_state=10)
        cluster_labels = clusterer.fit_predict(dataset)
        
        all_labels.append(cluster_labels)
        if n_clusters == 1:
            silhouette_avg = 0
        else:
            silhouette_avg = silhouette_score(dataset, cluster_labels)
        scores.append(silhouette_avg)
        distortions.append(clusterer.fit(dataset).inertia_)
        
    best = max_index(scores) + k_min
        
    for i in range(len(all_labels)):
        prefix = ''
        if (i + k_min == best):
            prefix = '_BEST'
            
        write_to_csv(dataset, all_labels[i], f'{args.output_path}/kmeans_with_' + str(i + k_min) + prefix + '_clusters.tsv')
        
        
        if (prefix == '_BEST'):
            labels = all_labels[i]
            predict = [x+1 for x in labels]
            classe = (pd.DataFrame(list(zip(dataset.index, predict)))).astype(str)
            classe.to_csv(best_cluster, sep = '\t', index = False, header = ['Patient_ID', 'Class'])
            
          
        
       
        if silhouette:
            silhouette_draw(dataset, all_labels[i], i + k_min, f'{args.output_path}/silhouette_with_' + str(i + k_min) + prefix + '_clusters.png')
        
        
    if elbow:
        elbow_plot(distortions, k_min,k_max) 

   
    
    

############################## elbow_plot ####################################
def elbow_plot (distortions: List[float], k_min: int, k_max: int) -> None:
    """
    Generate an elbow plot to visualize the distortion for different numbers of clusters.
    The elbow plot is a graphical tool used in clustering analysis to help identifying the appropriate number of clusters by looking for the point where the rate of decrease
    in distortion sharply decreases, indicating the optimal balance between model complexity and clustering quality.

    Args:
        distortions (list): List of distortion values for different numbers of clusters.
        k_min (int): The minimum number of clusters considered.
        k_max (int): The maximum number of clusters considered.

    Returns:
        None
    """
    plt.figure(0)
    x = list(range(k_min, k_max + 1))
    x.insert(0, 1)
    plt.plot(x, distortions, marker = 'o')
    plt.xlabel('Number of clusters (k)')
    plt.ylabel('Distortion')
    s = f'{args.output_path}/elbow_plot.png'
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward = True)
    fig.savefig(s, dpi=100)
    
    
############################## silhouette plot ###############################
def silhouette_draw(dataset: pd.DataFrame, labels: List[str], n_clusters: int, path:str) -> None:
    """
    Generate a silhouette plot for the clustering results.
    The silhouette coefficient is a measure used to evaluate the quality of clusters obtained from a clustering algorithmand it quantifies how similar an object is to its own cluster compared to other clusters.
    The silhouette coefficient ranges from -1 to 1, where:
    - A value close to +1 indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters. This implies that the object is in a dense, well-separated cluster.
    - A value close to 0 indicates that the object is close to the decision boundary between two neighboring clusters.
    - A value close to -1 indicates that the object may have been assigned to the wrong cluster.

    Args:
        dataset (pandas.DataFrame): The dataset used for clustering.
        labels (list): The cluster labels assigned to each data point.
        n_clusters (int): The number of clusters.
        path (str): The path to save the silhouette plot image.

    Returns:
        None
    """
    if n_clusters == 1:
        return None
        
    silhouette_avg = silhouette_score(dataset, labels)
    warning("For n_clusters = " + str(n_clusters) +
          " The average silhouette_score is: " + str(silhouette_avg))
           
    plt.close('all')
    # Create a subplot with 1 row and 2 columns
    fig, (ax1) = plt.subplots(1, 1)
    
    fig.set_size_inches(18, 7)
        
    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(dataset) + (n_clusters + 1) * 10])
    
    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(dataset, labels)
        
    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
        sample_silhouette_values[labels == i]
        
        ith_cluster_silhouette_values.sort()
    
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
    
        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                                     facecolor=color, edgecolor=color, alpha=0.7)
        
        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        
        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples
    
        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")
        
        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
    
        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
        
        
        plt.suptitle(("Silhouette analysis for clustering on sample data "
                          "with n_clusters = " + str(n_clusters) + "\nAverage silhouette_score = " + str(silhouette_avg)), fontsize=12, fontweight='bold')
            
            
        plt.savefig(path, bbox_inches='tight')
            
######################## dbscan ##############################################
def dbscan(dataset: pd.DataFrame, eps: float, min_samples: float, best_cluster: str) -> None:
    """
    Perform DBSCAN clustering on the given dataset, which is a clustering algorithm that groups together closely packed points based on the notion of density.

    Args:
        dataset (pandas.DataFrame): The dataset to be clustered.
        eps (float): The maximum distance between two samples for one to be considered as in the neighborhood of the other.
        min_samples (float): The number of samples in a neighborhood for a point to be considered as a core point.
        best_cluster (str): The file path to save the output of the best cluster.

    Returns:
        None
    """
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
        
    if eps is not None:
        clusterer = DBSCAN(eps = eps, min_samples = min_samples)
    else:
        clusterer = DBSCAN()
    
    clustering = clusterer.fit(dataset)
    
    core_samples_mask = np.zeros_like(clustering.labels_, dtype=bool)
    core_samples_mask[clustering.core_sample_indices_] = True
    labels = clustering.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    
    
    labels = labels
    predict = [x+1 for x in labels]
    classe = (pd.DataFrame(list(zip(dataset.index, predict)))).astype(str)
    classe.to_csv(best_cluster, sep = '\t', index = False, header = ['Patient_ID', 'Class'])
  
    
########################## hierachical #######################################
def hierachical_agglomerative(dataset: pd.DataFrame, k_min: int, k_max: int, best_cluster: str, silhouette: str) -> None:
    """
    Perform hierarchical agglomerative clustering on the given dataset.

    Args:
        dataset (pandas.DataFrame): The dataset to be clustered.
        k_min (int): The minimum number of clusters to consider.
        k_max (int): The maximum number of clusters to consider.
        best_cluster (str): The file path to save the output of the best cluster.
        silhouette (str): Whether to generate silhouette plots ('True' or 'False').

    Returns:
        None
    """
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    
    plt.figure(figsize=(10, 7))  
    plt.title("Customer Dendograms")  
    shc.dendrogram(shc.linkage(dataset, method='ward'), labels=dataset.index.values.tolist())  
    fig = plt.gcf()
    fig.savefig(f'{args.output_path}/dendogram.png', dpi=200)
    
    range_n_clusters = [i for i in range(k_min, k_max+1)]

    scores = []
    labels = []
    
    n_classi = dataset.shape[0]
    
    for n_clusters in range_n_clusters:  
        cluster = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage='ward')  
        cluster.fit_predict(dataset)  
        cluster_labels = cluster.labels_
        labels.append(cluster_labels)
        write_to_csv(dataset, cluster_labels, f'{args.output_path}/hierarchical_with_' + str(n_clusters) + '_clusters.tsv')
        
    best = max_index(scores) + k_min
    
    for i in range(len(labels)):
        prefix = ''
        if (i + k_min == best):
            prefix = '_BEST'
        if silhouette == 'true':
            silhouette_draw(dataset, labels[i], i + k_min, f'{args.output_path}/silhouette_with_' + str(i + k_min) + prefix + '_clusters.png')
     
    for i in range(len(labels)):
        if (i + k_min == best):
            labels = labels[i]
            predict = [x+1 for x in labels]
            classe = (pd.DataFrame(list(zip(dataset.index, predict)))).astype(str)
            classe.to_csv(best_cluster, sep = '\t', index = False, header = ['Patient_ID', 'Class'])
            
    
############################# main ###########################################
def main(args_in:List[str] = None) -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    """
    global args
    args = process_args(args_in)

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    
    #Data read
    
    X = read_dataset(args.input)
    X = X.iloc[:, 1:]
    X = pd.DataFrame.to_dict(X, orient='list')
    X = rewrite_input(X)
    X = pd.DataFrame.from_dict(X, orient = 'index')
    
    for i in X.columns:
        if any(val is None or np.isnan(val) for val in X[i]):
            X = X.drop(columns=[i])
            
    if args.k_max != None:
       numero_classi = X.shape[0]
       while args.k_max >= numero_classi:
          err = 'Skipping k = ' + str(args.k_max) + ' since it is >= number of classes of dataset'
          warning(err)
          args.k_max = args.k_max - 1
    
    
    if args.cluster_type == 'kmeans':
        kmeans(args.k_min, args.k_max, X, args.elbow, args.silhouette, args.best_cluster)
    
    if args.cluster_type == 'dbscan':
        dbscan(X, args.eps, args.min_samples, args.best_cluster)
        
    if args.cluster_type == 'hierarchy':
        hierachical_agglomerative(X, args.k_min, args.k_max, args.best_cluster, args.silhouette)
        
##############################################################################
if __name__ == "__main__":
    main()
