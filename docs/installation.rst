Installation
============

As of |today|, XYalign is currently under active development as we
prepare for our first major release.  Therefore, prerelease versions should be
considered unstable for the time being.

Obtaining XYalign
-----------------

In the meantime, if you're interested in using XYalign, you can obtain it by
cloning the Github repository::

	git clone https://github.com/WilsonSayresLab/XYalign/

Operating System
----------------

XYalign has been tested on a variety of Linux and Mac operating systems, but has
not been tested on Windows.  This isn't to say it won't work, however
we are unprepared to offer any Windows support at this time.

Requirements
------------

XYalign has a number of required Python packages and external programs::

	Python: 2.7

	Python packages:
		matplotlib
		numpy
		pandas
		pybedtools
		pysam
		seaborn

	External Programs:
		bbmap
		bedtools
		bwa
		platypus
		sambamba
		samtools

.. note::
	Bedtools is required for pybedtools and must be added to one's `PATH`. XYalign
	will check that it is available by calling `bedtools`. Other external programs
	do not, however, need to be on one's `PATH` and can be provided to XYalign
	using the appropriate flag(s)::

		--repairsh_path
		--bwa_path
		--platypus_path
		--samtools_path
		--sambamba_path

We strongly recommend users install and manage all packages and programs using
Anaconda.  To do so:

1. First download and install either
`Anaconda <https://www.continuum.io/downloads>`_
or `Miniconda <http://conda.pydata.org/miniconda.html>`_ (both work well,
Miniconda is a lightweight version of Anaconda).

	* Be sure to allow Anaconda to append to your PATH (it will ask for permission to do so during installation)

		* You can check this after installation with the command (from the command line)::

			which python

		which should point you to the python installed in your Anaconda or
		Miniconda directory.

2. Linux users can finish installation with the following commands (note that \\ indicates a continuation of the command on the next line)::

	conda config --add channels r

	conda config --add channels conda-forge

	conda config --add channels bioconda

	conda create -n xyalign_env python=2.7 pysam pybedtools \
	numpy pandas matplotlib seaborn platypus-variant bwa bbmap \
	samtools bedtools sambamba

And load the new environment (containing all required programs and packages) with::

	source activate xyalign_env

3. Mac users - as of right now, bioconda won't install platypus on Macs, so Mac
users will have to use the commands::

	conda config --add channels r

	conda config --add channels conda-forge

	conda config --add channels bioconda

	conda create -n xyalign_env python=2.7 pysam pybedtools \
	numpy pandas matplotlib seaborn bwa bbmap samtools bedtools sambamba

and then `install platypus on their own <http://www.well.ox.ac.uk/platypus>`_ and
provide it to XYalign with the flag::

	--platypus_path

Mac users can then load the environment with the command (same as Linux)::

		source activate xyalign_env
