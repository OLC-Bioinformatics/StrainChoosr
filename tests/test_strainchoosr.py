#!/usr/bin/env python

import pytest
import ete3
from strainchoosr.strainchoosr import *


def test_read_weights_file_good():
    weight_dict = read_weights_file('tests/text_files/good_weights_file.txt')
    assert weight_dict['strain1'] == 2.3
    assert weight_dict['strain2'] == 3.6
    assert len(weight_dict) == 2


def test_weights_file_wrong_separator():
    with pytest.raises(RuntimeError):
        read_weights_file('tests/text_files/weights_wrong_separator.txt')


def test_weights_file_bad_second_column():
    with pytest.raises(ValueError):
        read_weights_file('tests/text_files/weights_wrong_second_column.txt')


def test_weights_file_with_blank_lines():
    weight_dict = read_weights_file('tests/text_files/good_weights_file_with_blank_line.txt')
    assert weight_dict['strain1'] == 2.3
    assert weight_dict['strain2'] == 3.6
    assert len(weight_dict) == 2


def test_tree_modification_one_branch():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    weights = {'2018-SEQ-1315.fasta': 2}
    newtree = modify_tree_with_weights(tree, weights)
    # Original distance on that branch is 0.00002
    branch = newtree.get_leaves_by_name('2018-SEQ-1315.fasta')
    assert branch[0].dist == 0.00004
    assert len(branch) == 1
    # Also make sure original tree hasn't been modified
    branch = tree.get_leaves_by_name('2018-SEQ-1315.fasta')
    assert branch[0].dist == 0.00002


def test_tree_modification_multiple_branches():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    weights = {'2018-SEQ-1315.fasta': 2, '2018-SEQ-1271.fasta': 0.5}
    newtree = modify_tree_with_weights(tree, weights)
    # Original distance on that branch is 0.00002 - should be doubled.
    branch = newtree.get_leaves_by_name('2018-SEQ-1315.fasta')
    assert branch[0].dist == 0.00004
    # Original distance is 0.00003 - should be halved
    branch = newtree.get_leaves_by_name('2018-SEQ-1271.fasta')
    assert branch[0].dist == 0.000015


def test_tree_modification_bad_branch_name():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    weights = {'fake_branch': 3.6}
    with pytest.raises(AttributeError):
        modify_tree_with_weights(tree, weights)


def test_tree_modification_multiple_branches_same_name():
    tree = ete3.Tree('tests/tree_files/tree_multiple_same_name.nwk')
    weights = {'2018-SEQ-1315.fasta': 2}
    with pytest.raises(AttributeError):
        modify_tree_with_weights(tree, weights)


def test_starting_leaves_empty_list():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    starting_leaf_list = list()
    starting_leaves = find_starting_leaves(tree, starting_leaf_list)
    starter_names = list()
    for s in starting_leaves:
        starter_names.append(s.name)
    assert len(starter_names) == 2
    assert '2018-SEQ-0383.fasta' in starter_names
    assert '2018-SEQ-0100.fasta' in starter_names


def test_starting_leaves_one_starter():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    starting_leaf_list = list()
    starting_leaf = tree.get_leaves_by_name('2018-SEQ-0554.fasta')[0]
    starting_leaf_list.append(starting_leaf)
    starting_leaves = find_starting_leaves(tree, starting_leaf_list)
    starter_names = list()
    for s in starting_leaves:
        starter_names.append(s.name)
    assert len(starter_names) == 2
    assert '2018-SEQ-0383.fasta' in starter_names
    assert '2018-SEQ-0554.fasta' in starter_names


def test_starting_leaves_two_starters():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    starting_leaf_list = list()
    starting_leaf = tree.get_leaves_by_name('2018-SEQ-0554.fasta')[0]
    starting_leaf_list.append(starting_leaf)
    starting_leaf = tree.get_leaves_by_name('2018-SEQ-0525.fasta')[0]
    starting_leaf_list.append(starting_leaf)
    starting_leaves = find_starting_leaves(tree, starting_leaf_list)
    starter_names = list()
    for s in starting_leaves:
        starter_names.append(s.name)
    assert len(starter_names) == 2
    assert '2018-SEQ-0525.fasta' in starter_names
    assert '2018-SEQ-0554.fasta' in starter_names


def test_starting_leaves_more_than_two_starters():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    starting_leaf_list = list()
    starting_leaf = tree.get_leaves_by_name('2018-SEQ-0554.fasta')[0]
    starting_leaf_list.append(starting_leaf)
    starting_leaf = tree.get_leaves_by_name('2018-SEQ-0525.fasta')[0]
    starting_leaf_list.append(starting_leaf)
    starting_leaf = tree.get_leaves_by_name('2018-STH-0005.fasta')[0]
    starting_leaf_list.append(starting_leaf)
    starting_leaves = find_starting_leaves(tree, starting_leaf_list)
    starter_names = list()
    for s in starting_leaves:
        starter_names.append(s.name)
    assert len(starter_names) == 3
    assert '2018-SEQ-0525.fasta' in starter_names
    assert '2018-SEQ-0554.fasta' in starter_names
    assert '2018-STH-0005.fasta' in starter_names


def test_leaf_names_from_nodes():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    nodes = list()
    nodes.append(tree.get_leaves_by_name('2018-SEQ-0559.fasta')[0])
    nodes.append(tree.get_leaves_by_name('2018-SEQ-1315.fasta')[0])
    names = get_leaf_names_from_nodes(nodes)
    assert len(names) == 2
    assert '2018-SEQ-0559.fasta' in names
    assert '2018-SEQ-1315.fasta' in names


def test_find_next_leaf():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    starting_leaf_list = list()
    starting_leaves = find_starting_leaves(tree, starting_leaf_list)
    next_leaf = find_next_leaf(starting_leaves, tree)
    assert next_leaf.name == '2018-SEQ-0385.fasta'


def test_pd_greedy():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    starting_leaf_list = list()
    starting_leaves = find_starting_leaves(tree, starting_leaf_list)
    strains = pd_greedy('tests/tree_files/tree.nwk', 4, starting_leaves)
    assert len(strains) == 4
    names = get_leaf_names_from_nodes(strains)
    assert '2018-SEQ-0383.fasta' in names
    assert '2018-SEQ-0100.fasta' in names
    assert '2018-SEQ-0385.fasta' in names
    assert '2017-MER-0763.fasta' in names
