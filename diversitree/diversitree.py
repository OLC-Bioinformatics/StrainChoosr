#!/usr/bin/env python

import argparse
import logging
from Bio import Phylo
from scipy import cluster


def choose_representative_strain(cluster, tree):
    """
    Given a cluster, will find which tip in that cluster tends to be closest to all others (aka the best representative)
    and return that tip
    :param cluster: Cluster found.
    :param tree: Tree read in by Bio.Phylo
    :return: representative strain.
    """
    representative = 'NA'
    best_distance = 100000000.0  # Start off at an absolutely ridiculous value.
    # Iterate through
    for strain1 in cluster:
        total_length = 0.0
        for strain2 in cluster:
            if strain1 != strain2:
                total_length += tree.distance(strain1, strain2)
        if total_length < best_distance:
            best_distance = total_length
            representative = strain1
    return representative


def find_clusters(treefile, desired_clusters=10):
    """
    :param treefile: Tree file, in newick format.
    :param desired_clusters: The total number of clusters you want to return.
    :return: A list of strains that represent each cluster.
    """
    num_clusters = -1
    tree = Phylo.read(treefile, "newick")
    terminals = list()
    clades = tree.find_clades()
    matrix = list()
    # Get a list of all the terminal branches
    for clade in clades:
        if clade.is_terminal():
            terminals.append(str(clade))

    # Create a matrix with distances between each tip and each other tip
    for i in range(len(terminals)):
        for j in range(i + 1, len(terminals)):
            dist = tree.distance(terminals[i], terminals[j])
            matrix.append(dist)

    count = 1
    print('Finding correct number of clusters.')
    # Create the linkage thingy so we can try clustering.
    z = cluster.hierarchy.linkage(matrix, method='average')
    # Try different cutoff levels for making clusters. Start at a very large cluster distance.
    # Check if enough clusters are created. If not enough, make the cutoff slightly smaller and repeat process.
    cluster_distance = 0.9
    while num_clusters < desired_clusters:
        # clustering = cluster.hierarchy.fcluster(z, cluster_distance, criterion='distance')
        clustering = cluster.hierarchy.fcluster(z, cluster_distance, criterion='distance')
        num_clusters = (max(clustering))
        cluster_distance -= 0.00003
        count += 1
    # Once we've found the desired number of clusters, create a 2d array where each element of the array is a cluster
    # and each element is made up of a list of strains that are part of that cluster.
    clusters = list()
    for i in range(num_clusters):
        clusters.append(list())
    for i in range(len(clustering)):
        clusters[clustering[i] - 1].append(terminals[i])

    return clusters


def choose_representatives(clusters, tree):
    # Now that we have a list of strains for each cluster, choose the strain that's the most like everything else
    # in each cluster (the most average-y thing I guess?) to be the representative.
    strains = list()
    for c in clusters:
        representative_strain = choose_representative_strain(c, tree)
        representative_strain = representative_strain.replace('.fasta', '')
        strains.append(representative_strain)

    return strains


def draw_clustered_tree(treefile, clusters):
    colors = ['blue', 'red', 'green', 'purple', 'yellow', 'pink', 'gray']
    tree = Phylo.read(treefile, 'newick')
    i = 0
    for cluster in clusters:
        mrca = tree.common_ancestor(cluster)
        mrca.color = colors[i]
        i += 1
    Phylo.draw(tree)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--treefile',
                        type=str,
                        required=True,
                        help='Path to treefile, in newick format.')
    parser.add_argument('-n', '--number',
                        type=int,
                        required=True,
                        help='Number of representatives wanted.')
    args = parser.parse_args()
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    clusters = find_clusters(treefile=args.treefile,
                             desired_clusters=args.number)
    print(clusters)
    draw_clustered_tree(treefile=args.treefile,
                        clusters=clusters)
