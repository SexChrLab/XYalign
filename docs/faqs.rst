Frequently Asked Questions
==========================

Does XYalign require X and Y chromosomes?
---------------------------------------------

In principle, no, it doesn't.  The focus on X and Y chromosomes stems from our
initial interest in characterizing technical biases and aneuploidies affecting
variant calling on the sex chromosomes in large human genomic datasets. Hence,
the terminology we use throughout.  You can provide the name of any chromosome
or scaffold to ``--x_chromosome`` and ``--y_chromosome``, and an arbitrary number of
chromosome/scaffold names to ``--chromosomes``.  See :doc:`usage` for an example of
how this might work.  We plan to generalize XYalign in the future to make this
easier.

Will XYalign work with genomes from other organisms?
----------------------------------------------------

Yes, but with some caveats.  As discussed above, you can provide any chromosome
names to ``--x_chromosome`` and ``--y_chromosome``. So, if your organism
has Z and W chromosomes, this might look like ``--x_chromosome chrZ``
and ``--y_chromosome chrW``. However, we advise users to interpret results
cautiously, as XYalign's default settings for human X and Y chromosomes
are likely inappropriate for many other organisms.  This is especially the case
for ZW systems, or reference genomes without sequences for the Y (or equivalent)
chromosome.  In addition, XYalign does not currently accept multiple X or Y
scaffolds.  We plan to address these phenomena in future releases.
