|Travis| |Coverage| |PyPI| |ReadTheDocs|

StrainChoosr
============

``StrainChoosr`` examines phylogenetic trees and will give you the X strains that represent the most diversity
within your tree. Given today's deluge of sequencing data, picking strains to do more detailed analysis 
on can be difficult. Using this tool ensures you'll have the maximum amount of diversity possible in your
scaled down set of sequences. The algorithm behind it is described in Pardi2005_ and Steel2005_ so please be sure to cite them if you use
``StrainChoosr``.

Installation
============

``StrainChoosr`` lives on PyPi, so you can install via pip/pip3:

``pip install strainchoosr``

Alternatively, you can install the latest version from git. Note that there may be breaking changes
pushed to the repository, so this isn't necessarily advised.

``pip install git+https://github.com/lowandrew/StrainChoosr.git``

Quickstart
==========

To use ``StrainChoosr``, all you need is a newick-formatted_ tree file.
To pick the 5 most diverse strains from a tree, type the following on the command line:

``strainchoosr --treefile /path/to/tree.nwk --number 5``

This will print the names of the 5 most diverse strains in your tree to the terminal, as well as
create a file called `strainchoosr_output.html` in your current working directory that lets you visualize
the output in any web browser.

To do the same within python::

    from strainchoosr import strainchoosr
    strainchoosr.run_strainchoosr(treefile='/path/to/tree.nwk', number_representatives=[5])

In addition to printing the strains to terminal, ``run_strainchoosr`` will return a dictionary where keys are
the number of representatives and values are lists of the strains selected for that
number of representatives.

Alternatively, if all you want to get is the list of strains and not generated html reports::

    from strainchoosr import strainchoosr
    import ete3
    tree = ete3.Tree('path/to/treefile.nwk')
    diverse_strains = strainchoosr.pd_greedy(tree=tree, number_tips=5, starting_strains=[])

This will get you a list of ete3.TreeNode objects that represent the 5 most diverse possible strains. You can then use
``strainchoosr.get_leaf_names_from_nodes(diverse_strains)`` to get a list of names.

Complete documentation on the ``strainchoosr`` API can be found at https://strainchoosr.readthedocs.io/api.html.

Issues and Pull Requests
========================

If you have any problems or want a feature implemented, please feel free to open an issue_. Similarly, if you want to
add a feature or otherwise improve things, feel free to open a pull request.

.. _Pardi2005: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0010071
.. _Steel2005: https://academic.oup.com/sysbio/article/54/4/527/2842877
.. _newick-formatted: https://en.wikipedia.org/wiki/Newick_format)
.. _issue: https://github.com/lowandrew/StrainChoosr/issues
.. |Travis| image:: https://travis-ci.org/lowandrew/StrainChoosr.svg?master
           :target: https://travis-ci.org/lowandrew/StrainChoosr

.. |Coverage| image:: https://codecov.io/gh/lowandrew/StrainChoosr/branch/master/graph/badge.svg
           :target: https://codecov.io/gh/lowandrew/StrainChoosr

.. |PyPI| image:: https://badge.fury.io/py/strainchoosr.svg
           :target: https://badge.fury.io/py/strainchoosr
.. |ReadTheDocs| image:: https://readthedocs.org/projects/strainchoosr/badge/?version=latest
          :target: https://strainchoosr.readthedocs.io/?badge=latest
          :alt: Documentation Status