#!/usr/bin/env python

# Core python library
import os
import copy
import base64
import logging
import tempfile
import argparse

# Other stuff
import ete3
from ete3 import NodeStyle, TreeStyle, TextFace


def find_starting_leaves(tree, starting_leaf_list):
    """
    Gets the start of what the most diverse set of strains should be.
    If starting_leaf_list has nothing in it, the two strains that are the farthest apart will be picked (longest total
    branch length between them). If there is one entry in starting leaf list, the list returned will be that leaf plus
    whatever leaf on the tree is the farthest from it. If more than two leaves are in starting_leaf_list, nothing
    happens and just that will get returned.
    :param tree: An ete3.Tree object
    :param starting_leaf_list: A list of ete3.TreeNode objects.
    :return: A list of ete3.TreeNode objects representing the most diverse starting set possible
    """
    if len(starting_leaf_list) == 0:
        leaves = tree.get_leaves()
        max_distance = 0
        most_distant_leaves = None, None
        for leaf_one in leaves:
            for leaf_two in leaves:
                distance = tree.get_distance(leaf_one, leaf_two)
                if distance > max_distance:
                    max_distance = distance
                    most_distant_leaves = leaf_one, leaf_two
        starting_leaf_list.append(most_distant_leaves[0])
        starting_leaf_list.append(most_distant_leaves[1])
    elif len(starting_leaf_list) == 1:
        leaves = tree.get_leaves()
        max_distance = 0
        most_distant_leaf = None
        starting_leaf = starting_leaf_list[0]
        for leaf in leaves:
            if leaf.name != starting_leaf:
                distance = tree.get_distance(leaf, starting_leaf)
                if distance > max_distance:
                    most_distant_leaf = leaf
                    max_distance = distance
        starting_leaf_list.append(most_distant_leaf)

    return starting_leaf_list


def get_leaf_names_from_nodes(leaf_nodes):
    """

    :param leaf_nodes:
    :return:
    """
    leaf_names = list()
    for node in leaf_nodes:
        leaf_names.append(node.name)
    return leaf_names


def find_next_leaf(diverse_leaves, tree):
    """
    Given a set of leaves we've already decided represent the most diversity, find the next leaf that contributes
    the most diversity.
    :param diverse_leaves: List of leaves that we've already decided represent the most diversity possible - each
    entry in this list should be an ete3.TreeNode object
    :param tree: an ete3.Tree object that contains the nodes listed in diverse_leaves
    :return:
    """
    # Here, we prune off everything except for the leaves we've already selected as diverse and one other leaf in the tree
    # for each leaf in the tree, keeping branch lengths intact. Select the leaf that creates the tree with the longest
    # total branch length, and return that.
    longest_total_branch_length = 0
    sets_to_try = dict()
    leaves = tree.get_leaves()
    leaf_to_return = None
    for leaf in leaves:
        if leaf not in diverse_leaves:
            leafset = diverse_leaves.copy()
            leafset.append(leaf)
            sets_to_try[leaf.name] = leafset

    # In the event multiple strains have same distance, sort the keys so we're consistent about which one we're taking
    for leafname in sorted(sets_to_try):
        total_branch_length = 0
        newtree = tree.copy()
        newtree.prune(get_leaf_names_from_nodes(sets_to_try[leafname]), preserve_branch_length=True)
        for branch in newtree.get_descendants():
            total_branch_length += branch.dist
        if total_branch_length > longest_total_branch_length:
            leaf_to_return = tree.get_leaves_by_name(leafname)[0]
            longest_total_branch_length = total_branch_length
    return leaf_to_return


def pd_greedy(treefile, number_tips, starting_strains):
    # TODO here: implement greedy algorithm described in Species Choice for Comparative Genomics: Being Greedy Works
    #  (Pardi 2005), as well as the Steel 2005 paper
    # My understanding of this algorithm - start out by picking the two strains that have the longest total length
    # between them in the tree.
    # From there, add the leaf that adds the most total branch length to the tree, then just keep doing that until
    # you hit the number of strains you want.
    # To actually implement this, should just be able to have all strains but the ones in your tree pruned from the
    # original, and calculate total branch length after that.
    # Step 1: Read in ye olde treefile, figure out which nodes within tree structure are leaves.
    tree = ete3.Tree(newick=treefile)

    # Step 2: Find the two terminals with the most distance between them.
    starting_strains = find_starting_leaves(tree, starting_strains)

    while len(starting_strains) < number_tips:
        next_leaf = find_next_leaf(starting_strains, tree)
        if next_leaf is None:
            logging.error('ERROR: Could not find the next leaf to add.')
            quit(code=1)
        starting_strains.append(next_leaf)
    return starting_strains


def modify_tree_with_weights(tree, weights):
    """
    Given an ete3 Tree object and a dictionary where keys are node names in the tree and values are multipliers (can
    be generated with read_weights_file), returns a new tree where each branch in the weights dictionary is multiplied
    by the multiplier specified.
    :param tree: an ete3.Tree object
    :param weights: Dictionary where keys are names of nodes/tips in the tree, and values are weights by which branch
    lengths will be multiplied
    :return: A new ete3.Tree where branch lengths have been modified.
    """
    newtree = copy.deepcopy(tree)
    for node in weights:
        # Make sure that we can actually find the node, and that more than one branch doesn't have the same name.
        branch = newtree.get_leaves_by_name(node)
        if len(branch) != 1:
            raise AttributeError('The branch {} either could not be found in your tree or was found more than once. '
                                 'Please verify your tree/weights dictionary and try again.'.format(node))
        else:
            branch[0].dist *= weights[node]
    return newtree


def create_colored_tree_tip_image(tree_to_draw, representatives, output_file, color='red'):
    """
    Given a list of representatives, shows (for now) a phylogeny that has those representatives highlighted in
    a color to show it off.
    :param tree_to_draw: an ete3 Tree object. This won't get modified at any point.
    :param representatives: List with each strain name that should be highlighted.
    :param output_file: File to write output to, including extension.
    :param color: Color to show selected strains as. Defaults to red.
    """
    tree = copy.deepcopy(tree_to_draw)  # Don't want to actually modify original tree.
    ts = TreeStyle()
    ts.show_leaf_name = False
    for terminal_clade in tree.get_leaves():
        if terminal_clade.name in representatives:
            nstyle = NodeStyle()
            nstyle['shape'] = 'circle'
            nstyle['fgcolor'] = 'red'
            nstyle['size'] = 10
            name_face = TextFace(terminal_clade.name, fgcolor=color, fsize=10)
            terminal_clade.add_face(name_face, column=0)
            terminal_clade.set_style(nstyle)
        else:
            name_face = TextFace(terminal_clade.name, fgcolor='black', fsize=8)
            terminal_clade.add_face(name_face, column=0)

    tree.ladderize()
    tree.render(output_file, dpi=300, tree_style=ts)


def read_weights_file(weights_file):
    """
    Given a tab separated file with leaf names for a phylogenetic tree in column one and multipliers for that leaf's
    branch length in column two, will create a dictionary with leaf names as keys and multipliers as values
    :param weights_file: Path to a tab-separated text file described above.
    :return: dictionary with leaf names as keys and multipliers as values
    """
    weights = dict()
    with open(weights_file) as f:
        for line in f:
            stripped_line = line.rstrip()
            x = stripped_line.split('\t')
            print(x)
            if len(x) != 2 and stripped_line != '':
                raise RuntimeError('One of the lines in your weights file ({}) is not formatted correctly. '
                                   'Correct format is leafname\tweight, tab-separated. '
                                   'Offending line was: {}'.format(weights_file, stripped_line))
            elif len(x) == 2:
                try:
                    weight = float(x[1])
                except ValueError:
                    raise ValueError('The second column in your weights file ({}) must be a number. Please fix the '
                                     'following line: {}'.format(weights_file, stripped_line))
                weights[x[0]] = weight

    return weights


class CompletedStrainChoosr:
    def __init__(self, representatives, image, name):
        self.representatives = representatives
        self.image = image
        self.name = name


def generate_html_report(completed_choosr_list, output_report):
    # With tabs as shown in w3schools: https://www.w3schools.com/howto/howto_js_tabs.asp
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
    parser = argparse.ArgumentParser(description='StrainChoosr uses the greedy algorithm described in Pardi 2005/Steel '
                                                 '2005 to find the most diverse subset of strains from a phylogenetic '
                                                 'tree. ADD DOIs HERE')
    parser.add_argument('-t', '--treefile',
                        type=str,
                        required=True,
                        help='Path to treefile, in newick format.')
    parser.add_argument('-n', '--number',
                        type=int,
                        nargs='+',
                        required=True,
                        help='Number of representatives wanted. More than one can be specified, separated by '
                             'spaces.')
    parser.add_argument('-o', '--output_name',
                        default='strainchoosr_output',
                        type=str,
                        help='Base output name for file. PUT MORE INFO HERE.')
    parser.add_argument('--tree_mode',
                        default='r',
                        choices=['r', 'c'],
                        help='Mode to display output trees in - choose from r for rectangular ' 
                             'or c for circular. Defaults to rectangular.')
    parser.add_argument('--weight_file',
                        required=False,
                        help='Path to file specifying weights for leaves in tree. File must be tab-separated, with '
                             'leaf names in the first column and weights in the second. Leaves not listed will be '
                             'assigned a weight of 1.')
    parser.add_argument('--starting_strains',
                        default=list(),
                        nargs='+',
                        help='Names of strains that must be included in your set of diverse strains, separated by '
                             'spaces.')
    parser.add_argument('--verbosity',
                        choices=['debug', 'info', 'warning'],
                        default='info',
                        help='Choice of how much information you want printed to the terminal. Set debug to see a '
                             'ridiculous amount of stuff, info for a normal amount, and warning for very minimal '
                             'output.')
    args = parser.parse_args()

    if args.verbosity == 'info':
        logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                            level=logging.INFO,
                            datefmt='%Y-%m-%d %H:%M:%S')
    elif args.verbosity == 'debug':
        logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                            level=logging.DEBUG,
                            datefmt='%Y-%m-%d %H:%M:%S')
    elif args.verbosity == 'warning':
        logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                            level=logging.WARNING,
                            datefmt='%Y-%m-%d %H:%M:%S')
    # Now actually do stuff here!
    tree = ete3.Tree(args.treefile)
    starting_leaves = find_starting_leaves(tree, args.starting_strains)
    for number in args.number:
        strains = pd_greedy(args.treefile, number, starting_leaves)
        print(strains)


if __name__ == '__main__':
    main()
