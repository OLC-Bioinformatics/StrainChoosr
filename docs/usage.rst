Usage
=====

This section handles how to use ``StrainChoosr``. While it goes over command line options, the API
for the ``run_strainchoosr`` function is a mirror of the command line, so it should be applicable to both.

``StrainChoosr`` is built to find the most diverse subset of strains from a phylogenetic tree. To do this,
it uses a greedy algorithm described in Pardi2005_ and Steel2005_. Briefly, the algorithm starts by
finding the two leaves of a tree with the greatest distance between them, and then iteratively adds the leaf that adds
the most diversity possible to that tree until the desired number of strains has been found.

The basic usage is simple - from the command line the following will choose the 5 most diverse strains
from the newick-formatted treefile provided and print them to screen. Additionally, an HTML report (default name
`strainchoosr_output.html`) will be created that visualizes where on the tree your selected strains fall will be created.

``strainchoosr --treefile /path/to/tree.nwk --number 5``

It may be that you know for sure that you want certain strains selected, and go for the most diverse strains after that.
To do that, you can specify starting strains from the command line. For example, the following will start with `Strain1`
and `Strain3`, and add the 3 strains to those that create the most diversity.

``strainchoosr --treefile /path/to/tree.nwk --number 5 --starting_strains Strain1 Strain3``

Another possibility is that you have preferences for strains you want selected, but they aren't absolute - for this case,
you can make strains more or less likely to be chosen by modifying their weights. To do this, create a **tab-separated**
file with strain names in column one and the weight to modify by in column two - numbers greater than one increase
the chances a strain will be chosen, and numbers less than one decrease chances. This is done by multiplying branch length,
and so won't have any effect of 0 length branches.

The following weight file would make Strain1 more likely to be chosen and Strain3 less likely.::

    Strain1 2
    Strain3 0.2

You would then include this weight file (pretend you named it `weights.tsv`) in your comma would then include this weight file (pretend you named it `weights.tsv`) in your command line call.

``strainchoosr --treefile /path/to/tree.nwk --number 5 --weight_file weights.tsv``

A few other options that provide minor tweaks are available - full usage is below::

    usage: strainchoosr [-h] -t TREEFILE -n NUMBER [NUMBER ...] [-o OUTPUT_NAME]
                        [--tree_mode {r,c}] [--weight_file WEIGHT_FILE]
                        [--starting_strains STARTING_STRAINS [STARTING_STRAINS ...]]
                        [--color COLOR] [--verbosity {debug,info,warning}] [-v]

    StrainChoosr uses the greedy algorithm described in Pardi 2005/Steel 2005 to
    find the most diverse subset of strains from a phylogenetic tree.

    optional arguments:
      -h, --help            show this help message and exit
      -t TREEFILE, --treefile TREEFILE
                            Path to treefile, in newick format.
      -n NUMBER [NUMBER ...], --number NUMBER [NUMBER ...]
                            Number of representatives wanted. More than one can be
                            specified, separated by spaces.
      -o OUTPUT_NAME, --output_name OUTPUT_NAME
                            Base output name for file. PUT MORE INFO HERE.
      --tree_mode {r,c}     Mode to display output trees in - choose from r for
                            rectangular or c for circular. Defaults to
                            rectangular.
      --weight_file WEIGHT_FILE
                            Path to file specifying weights for leaves in tree.
                            File must be tab-separated, with leaf names in the
                            first column and weights in the second. Leaves not
                            listed will be assigned a weight of 1.
      --starting_strains STARTING_STRAINS [STARTING_STRAINS ...]
                            Names of strains that must be included in your set of
                            diverse strains, separated by spaces.
      --color COLOR         Color you want to have selected strains shown as. List
                            of available colors is available at http://etetoolkit.
                            org/docs/latest/reference/reference_treeview.html#ete3
                            .SVG_COLORS Defaults to red.
      --verbosity {debug,info,warning}
                            Choice of how much information you want printed to the
                            terminal. Set debug to see a ridiculous amount of
                            stuff, info for a normal amount, and warning for very
                            minimal output.
      -v, --version         show program's version number and exit


.. _Pardi2005: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0010071
.. _Steel2005: https://academic.oup.com/sysbio/article/54/4/527/2842877
