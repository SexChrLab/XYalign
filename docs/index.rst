XYalign: Inferring and Correcting for Sex Chromosome Ploidy in Next-Generation Sequencing Data
==============================================================================================
:Authors: Timothy Webster, Tanya Phung, Madeline Couse, Bruno Grande, Eric Karlins, Phillip Richmond, Whitney Whitford, Melissa Wilson Sayres
:Date: |today|
:Release: |release|
:Download: `Github Repository <https://github.com/WilsonSayresLab/XYalign/>`_

Sex chromosome aneuploidies are currently estimated to be as common as 1/400 in humans. Atypical ploidy will affect variant calling and measures of genomic variation that are central to most clinical genomic studies. Further, the high degree of similarity between gametologous sequences on the X and Y chromosomes can lead to the misalignment of sequencing reads and substantially affect variant calling. Here we present XYalign, a new tool that (1) quickly infers sex chromosome ploidy in NGS data (DNA and RNA), (2) remaps reads based on the inferred sex chromosome complement of the individual, and (3) outputs quality, depth, and allele-balance metrics across chromosomes.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation.rst
   usage.rst
   faqs.rst
   API </api/modules>
   release.rst




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
