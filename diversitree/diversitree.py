#!/usr/bin/env python

import argparse
import logging
from Bio import Phylo
from scipy import cluster
import matplotlib.pyplot as plt


class DiversiTree(object):
    def __init__(self, tree_file, tree_format='newick'):
        self.tree = Phylo.read(tree_file, tree_format)
        self.terminal_clades = self.tree.get_terminals()

    def create_linkage(self, linkage_method='average'):
        matrix = list()
        for i in range(len(self.terminal_clades)):
            for j in range(i + 1, len(self.terminal_clades)):
                distance_between_tips = self.tree.distance(self.terminal_clades[i], self.terminal_clades[j])
                matrix.append(distance_between_tips)

        return cluster.hierarchy.linkage(matrix, method=linkage_method)

    def draw_colored_dendrogram(self, linkage):
        # TODO: Figure this out more.
        fig, axes = plt.subplots(1, 1)
        dn = cluster.hierarchy.dendrogram(linkage, ax=axes)
        plt.show()

    def find_clusters(self, linkage, desired_clusters=5, criterion='distance', step_size=0.00003):
        # Initialize our values - set initial cluster distance to ridiculous value, and num clusters too low.
        cluster_distance = 0.9
        num_clusters = -1
        while num_clusters < desired_clusters:
            clustering = cluster.hierarchy.fcluster(linkage, cluster_distance, criterion=criterion)
            num_clusters = max(clustering)
            cluster_distance -= step_size

        # Once we have desired number of clusters, create a 2d array where each element of the array is a cluster
        # and each element is made up of a list of strains that are part of that cluster.
        clusters = list()
        for i in range(num_clusters):
            clusters.append(list())
        for i in range(len(clustering)):
            clusters[clustering[i] - 1].append(self.terminal_clades[i])
        return clusters

    def choose_best_representative(self, cluster, method='closest'):
        representative = 'NA'
        assert method in ('closest', 'farthest')
        if method == 'closest':
            best_distance = 100000000.0  # Start off at an absolutely ridiculous value.
        elif method == 'farthest':
            best_distance = -0.1
        # Iterate through
        for strain1 in cluster:
            total_length = 0.0
            for strain2 in cluster:
                if strain1 != strain2:
                    total_length += self.tree.distance(strain1, strain2)
            if method == 'closest':
                if total_length < best_distance:
                    best_distance = total_length
                    representative = strain1
            elif method == 'farthest':
                if total_length > best_distance:
                    best_distance = total_length
                    representative = strain1
        return representative

    def cluster_contains_root(self, cluster):
        mrca = self.tree.common_ancestor(cluster)
        if mrca == self.tree.root:
            return True
        else:
            return False


def draw_clustered_tree(treefile, clusters):
    # TODO: Needs extremely major renovations - this should save to a file, not show up and halt the script.
    # Also, need to get a better color scheme created.
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
    diversitree = DiversiTree(tree_file=args.treefile)
    linkage = diversitree.create_linkage()
    clusters = diversitree.find_clusters(linkage=linkage, desired_clusters=args.number)
    for cluster in clusters:
        print(diversitree.choose_best_representative(cluster, method='closest'))
        # print(diversitree.cluster_contains_root(cluster))
    # clusters = find_clusters(treefile=args.treefile,
    #                          desired_clusters=args.number)
    # for i in range(len(clusters)):
    #     print('Cluster {}'.format(i))
    #     print(choose_representative_strain(clusters[i], Phylo.read(args.treefile, 'newick')))
    # draw_clustered_tree(treefile=args.treefile,
    #                     clusters=clusters)
