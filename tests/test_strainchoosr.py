#!/usr/bin/env python

import pytest
import tempfile
import ete3
import os
from unittest.mock import patch
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


def test_leaf_nodes_from_names():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    tree_nodes = get_leaf_nodes_from_names(tree, ['2018-SEQ-0559.fasta', '2018-SEQ-1315.fasta'])
    assert len(tree_nodes) == 2


def test_leaf_nodes_from_names_bad_name():
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    with pytest.raises(RuntimeError):
        tree_nodes = get_leaf_nodes_from_names(tree, ['2018-SEQ-0559.fasta', 'super_fake_leaf'])


def test_get_version():
    # Todo: this test stinks. Make it more useful
    version = get_version()
    assert 'StrainChoosr' in version


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
    strains = pd_greedy(tree, 4, starting_leaves)
    assert len(strains) == 4
    names = get_leaf_names_from_nodes(strains)
    assert '2018-SEQ-0383.fasta' in names
    assert '2018-SEQ-0100.fasta' in names
    assert '2018-SEQ-0385.fasta' in names
    assert '2017-MER-0763.fasta' in names
    assert len(starting_leaves) == 2


def test_tree_draw():
    with tempfile.TemporaryDirectory() as tmpdir:
        tree = ete3.Tree('tests/tree_files/tree.nwk')
        representatives = ['2018-SEQ-0383.fasta', '2018-SEQ-0100.fasta', '2018-SEQ-0385.fasta', '2017-MER-0763.fasta']
        output_file = os.path.join(tmpdir, 'tree.png')
        create_colored_tree_tip_image(tree, representatives, output_file)
        assert os.path.isfile(output_file)


def test_tree_draw_pdf():
    with tempfile.TemporaryDirectory() as tmpdir:
        tree = ete3.Tree('tests/tree_files/tree.nwk')
        representatives = ['2018-SEQ-0383.fasta', '2018-SEQ-0100.fasta', '2018-SEQ-0385.fasta', '2017-MER-0763.fasta']
        output_file = os.path.join(tmpdir, 'tree.pdf')
        create_colored_tree_tip_image(tree, representatives, output_file)
        assert os.path.isfile(output_file)


def test_tree_draw_alternate_color():
    with tempfile.TemporaryDirectory() as tmpdir:
        tree = ete3.Tree('tests/tree_files/tree.nwk')
        representatives = ['2018-SEQ-0383.fasta', '2018-SEQ-0100.fasta', '2018-SEQ-0385.fasta', '2017-MER-0763.fasta']
        output_file = os.path.join(tmpdir, 'tree.png')
        create_colored_tree_tip_image(tree, representatives, output_file, color='lavender')
        assert os.path.isfile(output_file)


def test_tree_draw_rotated():
    with tempfile.TemporaryDirectory() as tmpdir:
        tree = ete3.Tree('tests/tree_files/tree.nwk')
        representatives = ['2018-SEQ-0383.fasta', '2018-SEQ-0100.fasta', '2018-SEQ-0385.fasta', '2017-MER-0763.fasta']
        output_file = os.path.join(tmpdir, 'tree.png')
        create_colored_tree_tip_image(tree, representatives, output_file, color='lavender', rotation=90)
        assert os.path.isfile(output_file)


def test_tree_draw_alternate_shape():
    with tempfile.TemporaryDirectory() as tmpdir:
        tree = ete3.Tree('tests/tree_files/tree.nwk')
        representatives = ['2018-SEQ-0383.fasta', '2018-SEQ-0100.fasta', '2018-SEQ-0385.fasta', '2017-MER-0763.fasta']
        output_file = os.path.join(tmpdir, 'tree.png')
        create_colored_tree_tip_image(tree, representatives, output_file, mode='c')
        assert os.path.isfile(output_file)


def test_html_report_generation():
    completed_choosrs = list()
    number = 4
    tree = ete3.Tree('tests/tree_files/tree.nwk')
    with tempfile.TemporaryDirectory() as tmpdir:
        representatives = ['2018-SEQ-0383.fasta', '2018-SEQ-0100.fasta', '2018-SEQ-0385.fasta', '2017-MER-0763.fasta']
        output_image = os.path.join(tmpdir, 'strains_{}.png'.format(number))
        create_colored_tree_tip_image(tree_to_draw=tree,
                                      output_file=output_image,
                                      representatives=representatives)
        completed_choosrs.append(CompletedStrainChoosr(representatives,
                                                       image=output_image,
                                                       name='{} Strains'.format(number)))
        generate_html_report(completed_choosrs,
                             os.path.join(tmpdir, 'strainchoosr_report.html'))


def test_argument_parsing_mostly_defaults():
    args = argument_parsing(['-t', 'tests/tree_files/tree.nwk', '-n', '5', '10', '20'])
    assert args.treefile == 'tests/tree_files/tree.nwk'
    assert args.number == [5, 10, 20]
    assert args.output_name == 'strainchoosr_output'
    assert args.tree_mode == 'r'
    assert args.weight_file is None
    assert args.starting_strains == []
    assert args.verbosity == 'info'


def test_argument_parsing_starting_strains():
    args = argument_parsing(['-t', 'tests/tree_files/tree.nwk', '-n', '5', '10', '20',
                             '--starting_strains', '2018-SEQ-0100.fasta'])
    assert args.treefile == 'tests/tree_files/tree.nwk'
    assert args.number == [5, 10, 20]
    assert args.output_name == 'strainchoosr_output'
    assert args.tree_mode == 'r'
    assert args.weight_file is None
    assert args.starting_strains == ['2018-SEQ-0100.fasta']
    assert args.verbosity == 'info'


def test_run_strainchoosr():
    with tempfile.TemporaryDirectory() as tmpdir:
        run_strainchoosr(treefile='tests/tree_files/tree.nwk',
                         number_representatives=[5, 10],
                         output_name=os.path.join(tmpdir, 'st_report'))
        assert os.path.isfile(os.path.join(tmpdir, 'st_report.html'))


def test_run_strainchoosr_output_dictionary():
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dict = run_strainchoosr(treefile='tests/tree_files/tree.nwk',
                                       number_representatives=[4],
                                       output_name=os.path.join(tmpdir, 'st_report'))
        assert os.path.isfile(os.path.join(tmpdir, 'st_report.html'))
    assert output_dict[4] == ['2018-SEQ-0383.fasta', '2018-SEQ-0100.fasta', '2018-SEQ-0385.fasta', '2017-MER-0763.fasta']


def test_run_strainchoosr_too_many_strains():
    with pytest.raises(ValueError):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_strainchoosr(treefile='tests/tree_files/tree.nwk',
                             number_representatives=[5, 10, 256],
                             output_name=os.path.join(tmpdir, 'st_report'))
            assert os.path.isfile(os.path.join(tmpdir, 'st_report.html'))


def test_run_strainchoosr_weights_file():
    with tempfile.TemporaryDirectory() as tmpdir:
        run_strainchoosr(treefile='tests/tree_files/tree.nwk',
                         number_representatives=[5, 10],
                         output_name=os.path.join(tmpdir, 'st_report'),
                         weight_file='tests/text_files/weights.txt')
        assert os.path.isfile(os.path.join(tmpdir, 'st_report.html'))


def test_run_strainchoosr_debug_log():
    with tempfile.TemporaryDirectory() as tmpdir:
        run_strainchoosr(treefile='tests/tree_files/tree.nwk',
                         number_representatives=[5, 10],
                         output_name=os.path.join(tmpdir, 'st_report'),
                         verbosity='debug')
        assert os.path.isfile(os.path.join(tmpdir, 'st_report.html'))


def test_run_strainchoosr_warning_log():
    with tempfile.TemporaryDirectory() as tmpdir:
        run_strainchoosr(treefile='tests/tree_files/tree.nwk',
                         number_representatives=[5, 10],
                         output_name=os.path.join(tmpdir, 'st_report'),
                         verbosity='warning')
        assert os.path.isfile(os.path.join(tmpdir, 'st_report.html'))


def test_completed_choosr_object():
    asdf = CompletedStrainChoosr(name='asdf',
                                 image='asdf.png',
                                 representatives=['asdf', 'fdsa'])
    assert asdf.name == 'asdf'
    assert asdf.image == 'asdf.png'
    assert asdf.representatives == ['asdf', 'fdsa']


def test_main():
    with tempfile.TemporaryDirectory() as tmpdir:
        output_stuff = os.path.join(tmpdir, 'st_output')
        fake_args = ['strainchoosr', '-t', 'tests/tree_files/tree.nwk', '-n', '5', '-o', output_stuff]
        with patch('sys.argv', fake_args):
            main()
            assert os.path.isfile(output_stuff + '.html')