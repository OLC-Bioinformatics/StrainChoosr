#!/usr/bin/env python

import glob
import shutil
import logging


def choose_tree_maker(fasta_folder, fasta_extension='.fasta'):
    """
    Given a folder of fasta-formatted files, try to do some auto figuring out regarding what program should
    be used to create the tree.
    Options:
    - If fasta files look to be gene sequences, do alignment (mafft/muscle/whatever) and then make a tree
    from that alignment (fasttree?)
    - If we have what appear to be a bunch of whole genome sequences, run mash to try to find out how far apart they
    tend to be. If they're close, use parsnp, if farther, mashtree.
    TODO: Handle FASTQ files? Compressed FASTA/FASTQ/other?

    :param fasta_folder:
    :param fasta_extension:
    :return:
    """
    return True

