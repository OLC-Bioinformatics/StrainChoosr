#!/usr/bin/env python

# Core python library
import os
import copy
import base64
import random
import logging
import tempfile
import argparse

# Other stuff
import ete3
from ete3 import NodeStyle, TreeStyle, TextFace
from scipy.cluster import hierarchy


def random_color():
    return '#' + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])


class StrainChoosr(object):
    def __init__(self, tree_file):
        self.tree = ete3.Tree(newick=tree_file)   # Would have as PhyloTree, but that makes causes problems with tree drawing
        self.terminal_clades = self.tree.get_leaves()

    def create_linkage(self, linkage_method='average'):
        """
        Finds the distance between each tip in a tree and puts it all into a matrix. Once the matrix has been created,
        does some hierarchical clustering.
        :param linkage_method: Method you want to use to perform clustering. Options are single, complete, average,
        weighted, centroid, median and ward. Defaults to 'average'. Lots more information on this available at
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
        :return: the result of scipy.cluster.hierarchy.linkage
        """
        linkage_methods = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']
        assert linkage_method in linkage_methods
        matrix = list()
        for i in range(len(self.terminal_clades)):
            for j in range(i + 1, len(self.terminal_clades)):
                distance_between_tips = self.tree.get_distance(target=self.terminal_clades[i], target2=self.terminal_clades[j])
                matrix.append(distance_between_tips)

        return hierarchy.linkage(matrix, method=linkage_method)

    def find_clusters(self, linkage, desired_clusters=5, criterion='distance', step_size=0.00003):
        """
        Creates clusters based on a linkage matrix. First creates with a very large threshold, decreasing threshold
        until the desired number of clusters is found.
        :param linkage: Linkage matrix created by create_linkage
        :param desired_clusters: Number of clusters you want your data to have.
        :param criterion: The criterion for forming clusters form linkage. Defaults to 'distance', which is the only
        one I've really tried. Other options are inconsistent, maxclust, monocrit, and maxclust_monocrit
        :param step_size: How much to decrease the cluster distance by each iteration. Too low and this takes a long
        time, too high and the number of desired clusters will probably get overshot.
        :return: A 2D array where each cluster is an array, and the names of the terminal clades are stored in array.
        For example, with two clusters with A, B, and C in cluster one and D and E in cluster 2, will return:
        [[A, B, C], [D, E]]
        """
        # Initialize our values - set initial cluster distance to ridiculous value, and num clusters too low.
        cluster_distance = 0.9
        num_clusters = -1
        while num_clusters < desired_clusters:
            clustering = hierarchy.fcluster(linkage, cluster_distance, criterion=criterion)
            num_clusters = max(clustering)
            cluster_distance -= step_size

        # Once we have desired number of clusters, create a 2d array where each element of the array is a cluster
        # and each element is made up of a list of strains that are part of that cluster.
        clusters = list()
        for i in range(num_clusters):
            clusters.append(list())
        for i in range(len(clustering)):
            clusters[clustering[i] - 1].append(self.terminal_clades[i].name)
        return clusters

    def choose_best_representative(self, cluster, method='closest'):
        """
        Given a cluster, will find the tip that best represents that cluster by either choosing the tip that is on
        average closest to the other tips in that cluster (the default) or whatever is on average farthest from other
        tips
        :param cluster: A list of terminal node names in the tree used to instantiate the object. Can be one of the
        clusters generated by find_clusters
        :param method: How to choose best representative. Either finds the tip that is on average closest to or farthest
        from all other tips in the cluster
        :return: Name of the tip chosen as best
        """
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
                    total_length += self.tree.get_distance(target=strain1, target2=strain2)
            if method == 'closest':
                if total_length < best_distance:
                    best_distance = total_length
                    representative = strain1
            elif method == 'farthest':
                if total_length > best_distance:
                    best_distance = total_length
                    representative = strain1
        return representative

    def create_colored_tree_tip_image(self, representatives, output_file):
        """
        Given a list of representatives, shows (for now) a phylogeny that has those representatives highlighted in
        a color to show it off.
        :param representatives: List with each strain name that should be highlighted.
        :param output_file: File to write output to, including extension.
        """
        tree = copy.deepcopy(self.tree)  # Don't want to actually modify original tree.
        ts = TreeStyle()
        ts.show_leaf_name = False
        for terminal_clade in tree.get_leaves():
            if terminal_clade.name in representatives:
                nstyle = NodeStyle()
                nstyle['shape'] = 'circle'
                nstyle['fgcolor'] = 'red'
                nstyle['size'] = 10
                name_face = TextFace(terminal_clade.name, fgcolor='red', fsize=10)
                terminal_clade.add_face(name_face, column=0)
                terminal_clade.set_style(nstyle)
            else:
                name_face = TextFace(terminal_clade.name, fgcolor='black', fsize=8)
                terminal_clade.add_face(name_face, column=0)

        tree.ladderize()
        tree.render(output_file, dpi=300, tree_style=ts)

    @staticmethod
    def find_common_ancestor_and_descendants(tree, cluster):
        nodes = list()
        for terminal in tree.get_leaves():
            if terminal.name in cluster:
                nodes.append(terminal)
        if len(nodes) == 1:
            return nodes[0]
        common_ancestor = nodes[0].get_common_ancestor(nodes[1:])
        descendants = common_ancestor.get_descendants()
        descendants.append(common_ancestor)  # Have to include this so that ancestor actually ends up colored.
        descendants.append(nodes[0])  # This doesn't get included in the get_descendants? I may have misunderstood docs
        # Check that the root of the tree isn't actually in the list - if it is, the whole tree ends up colored
        # one color, which is not at all what we want.
        root = tree.get_tree_root()
        if root in descendants:
            descendants.remove(root)
        return descendants

    def draw_clustered_tree(self, clusters, output_file):
        tree = copy.deepcopy(self.tree)
        common_ancestor_groups = list()
        for cluster in clusters:
            common_ancestor_groups.append(self.find_common_ancestor_and_descendants(tree, cluster))
        for common_ancestor_group in common_ancestor_groups:
            group_color = random_color()
            nstyle = NodeStyle()
            nstyle['hz_line_color'] = group_color
            nstyle['vt_line_color'] = group_color
            nstyle['fgcolor'] = group_color
            for node in common_ancestor_group:
                node.set_style(nstyle)
        tree.ladderize()
        if output_file is not None:
            tree.render(output_file, dpi=300)
        return tree

    def draw_clusters_and_tips(self, clusters, output_file, representatives):
        tree = self.draw_clustered_tree(clusters=clusters, output_file=None)
        ts = TreeStyle()
        ts.show_leaf_name = False
        for terminal_clade in tree.get_leaves():
            if terminal_clade.name in representatives:
                nstyle = NodeStyle()
                nstyle['shape'] = 'circle'
                nstyle['fgcolor'] = 'red'
                nstyle['size'] = 10
                name_face = TextFace(terminal_clade.name, fgcolor='red', fsize=10)
                terminal_clade.add_face(name_face, column=0)
                terminal_clade.set_style(nstyle)
            else:
                name_face = TextFace(terminal_clade.name, fgcolor='black', fsize=8)
                terminal_clade.add_face(name_face, column=0)
        tree.render(output_file, dpi=300, tree_style=ts)


class CompletedStrainChoosr:
    def __init__(self, representatives, image, name):
        self.representatives = representatives
        self.image = image
        self.name = name


def generate_html_report(completed_choosr_list, output_report):
    style = """
    <style>
    body {font-family: Arial;}

    /* Style the tab */
    .tab {
      overflow: hidden;
      border: 1px solid #ccc;
      background-color: #f1f1f1;
    }

    /* Style the buttons inside the tab */
    .tab button {
      background-color: inherit;
      float: left;
      border: none;
      outline: none;
      cursor: pointer;
      padding: 14px 16px;
      transition: 0.3s;
      font-size: 17px;
    }

    /* Change background color of buttons on hover */
    .tab button:hover {
      background-color: #ddd;
    }

    /* Create an active/current tablink class */
    .tab button.active {
      background-color: #ccc;
    }

    /* Style the tab content */
    .tabcontent {
      display: none;
      padding: 6px 12px;
      border: 1px solid #ccc;
      border-top: none;
    }
    </style>
    """

    javascript = """
    <script>
    function openCity(evt, cityName) {
      var i, tabcontent, tablinks;
      tabcontent = document.getElementsByClassName("tabcontent");
      for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = "none";
      }
      tablinks = document.getElementsByClassName("tablinks");
      for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" active", "");
      }
      document.getElementById(cityName).style.display = "block";
      evt.currentTarget.className += " active";
    }
    </script>
    """
    html_content = list()
    html_content.append('<html><head>')
    html_content.append(style)
    html_content.append('</head><body>')
    html_content.append('<h1>StrainChoosr Report</h1><br>')

    html_content.append('<div class="tab">\n')
    for completed_choosr in completed_choosr_list:
        html_content.append('<button class="tablinks" onclick="openCity(event, \'{name}\')">{name}</button>'.format(name=completed_choosr.name))
    html_content.append('</div>')
    for completed_choosr in completed_choosr_list:
        html_content.append('<div id="{name}" class="tabcontent">'.format(name=completed_choosr.name))
        html_content.append('<h4>{}</h4><br>'.format(completed_choosr.name))
        with open(completed_choosr.image, 'rb') as image_file:
            base64_string = base64.b64encode(image_file.read()).decode('utf-8')
        html_content.append('<img src="data:image/png;base64,{}">'.format(base64_string))
        html_content.append('<br><h4>Chosen Strains</h4>')
        for strain in completed_choosr.representatives:
            html_content.append('<p>{}</p>'.format(strain))
        html_content.append('</div>')

    html_content.append(javascript)
    html_content.append('</body></html>')
    html_string = '\n'.join(html_content)
    with open(output_report, 'w') as f:
        f.write(html_string)


def main():
    parser = argparse.ArgumentParser(description='StrainChoosr provides a set of tools for choosing the most diverse '
                                                 'set of strains from a set of DNA sequences or protein sequences.')
    subparsers = parser.add_subparsers(help='asdf', dest='subparsers')

    # SUBPARSER FOR TREE GENERATION
    treefile_subparser = subparsers.add_parser('choose', help='For strain choosing if you already have a tree file.')
    treefile_subparser.add_argument('-t', '--treefile',
                                    type=str,
                                    required=True,
                                    help='Path to treefile, in newick format.')
    treefile_subparser.add_argument('-n', '--number',
                                    type=int,
                                    nargs='+',
                                    required=True,
                                    help='Number of representatives wanted.')
    treefile_subparser.add_argument('-o', '--output_folder',
                                    type=str,
                                    required=True,
                                    help='Output folder to store results of StrainChoosr.')

    # SUBPARSER FOR TREE CREATION
    tree_creation = subparsers.add_parser('tree_create', help='Create a tree from a set of FASTA sequences.')
    tree_creation.add_argument('arg1')

    #$ SUBPARSER FOR CREATION AND CHOOSING
    create_and_choose = subparsers.add_parser('create_and_choose', help='Both creates a tree and chooses strains '
                                                                        'from the tree created.')
    create_and_choose.add_argument('asdf')

    # Actually start doing things.
    args = parser.parse_args()
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    # JUST CHOOSE STRAINS
    if args.subparsers == 'choose':
        if not os.path.isdir(args.output_folder):
            os.makedirs(args.output_folder)
        with tempfile.TemporaryDirectory() as tmpdir:
            completed_choosrs = list()
            diversitree = StrainChoosr(tree_file=args.treefile)
            linkage = diversitree.create_linkage()
            for number in args.number:
                clusters = diversitree.find_clusters(linkage=linkage, desired_clusters=number)
                reps = list()
                for cluster in clusters:
                    rep = diversitree.choose_best_representative(cluster, method='closest')
                    reps.append(rep)
                output_file = os.path.join(tmpdir, 'strains_{}.png'.format(number))
                diversitree.draw_clusters_and_tips(clusters=clusters, output_file=output_file, representatives=reps)
                completed_choosrs.append(CompletedStrainChoosr(representatives=reps,
                                                               image=output_file,
                                                               name='{} Strains'.format(number)))
            html_report = os.path.join(args.output_folder, 'StrainChoosr.html')
            generate_html_report(completed_choosr_list=completed_choosrs,
                                 output_report=html_report)

    # JUST CREATE A TREE
    elif args.subparsers == 'tree_create':
        print('Create tree ')

    # CREATE A TREE AND CHOOSE STRAINS


if __name__ == '__main__':
    main()

