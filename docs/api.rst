API
===

Reference/Fasta files
---------------------

.. autoclass:: xyalign.reftools.RefFasta
	:members:

Bam files
---------

.. autoclass:: xyalign.bam.BamFile
	:members:

.. autofunction:: xyalign.bam.sambamba_merge

.. autofunction:: xyalign.bam.switch_sex_chromosomes_sambamba

Variant Calling/ VCF Processing
-------------------------------

.. autofunction:: xyalign.variants.platypus_caller

.. autofunction:: xyalign.variants.parse_platypus_VCF

.. autofunction:: xyalign.variants.plot_variants_per_chrom

.. autofunction:: xyalign.variants.hist_read_balance

.. autofunction:: xyalign.variants.plot_read_balance

Ploidy Estimation
-----------------

.. autofunction:: xyalign.ploidy.bootstrap

.. autofunction:: xyalign.ploidy.ks_two_sample

.. autofunction:: xyalign.ploidy.permutation_test_chromosomes

Mapping/assembly
----------------

.. autofunction:: xyalign.assemble.bwa_mem_mapping_sambamba

Utilities
---------

.. autofunction:: xyalign.utils.check_bam_fasta_compatibility

.. autofunction:: xyalign.utils.chromosome_wide_plot

.. autofunction:: xyalign.utils.chromosome_bed

.. autofunction:: xyalign.utils.make_region_lists

.. autofunction:: xyalign.utils.merge_bed_files

.. autofunction:: xyalign.utils.output_bed

.. autofunction:: xyalign.utils.plot_depth_mapq
