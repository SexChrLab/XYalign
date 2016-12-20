API
===

Reference/Fasta files
---------------------

.. autoclass:: reftools.RefFasta
	:members:

Bam files
---------

.. autoclass:: bam.BamFile
	:members:

.. autofunction:: bam.sambamba_merge

.. autofunction:: bam.switch_sex_chromosomes_sambamba

Variant Calling/ VCF Processing
-------------------------------

.. autofunction:: variants.platypus_caller

.. autofunction:: variants.parse_platypus_VCF

.. autofunction:: variants.plot_variants_per_chrom

.. autofunction:: variants.hist_read_balance

.. autofunction:: variants.plot_read_balance

Ploidy Estimation
-----------------

.. autofunction:: ploidy.bootstrap

.. autofunction:: ploidy.ks_two_sample

.. autofunction:: ploidy.permutation_test_chromosomes

Mapping/assembly
----------------

.. autofunction:: assemble.bwa_mem_mapping_sambamba

Utilities
---------

.. autofunction:: utils.check_bam_fasta_compatibility

.. autofunction:: utils.chromosome_wide_plot

.. autofunction:: utils.chromosome_bed

.. autofunction:: utils.make_region_lists

.. autofunction:: utils.merge_bed_files

.. autofunction:: utils.output_bed

.. autofunction:: utils.plot_depth_mapq
