A snakemake pipeline for identification of extreme reads pile-up sites across genome from 10x Genomics single cell RNA-Seq experiment.

Installation
------------

* Step1: install miniconda3 of your system
.. code:: bash

   $ sh Miniconda3-latest-MacOSX-x86_64.sh


* Step2: install snakemake
.. code:: bash

   $ conda install -c conda-forge -c bioconda snakemake=5.11.2



   Check snakemake installed successfully

   $ which snakemake && snakemake --version

   if you didn't see snakemake in your path, you may close and re-open your terminal to activate the effects, and check again.

 * Step3: download and run
 .. code:: bash

    $ git clone https://github.com/yh154/Pyup.github
    $ cd Pyup
    $ snakemake --cores 1 --use-conda -s Snakefile

Input and Output
----------------
   The pipeline expects a single bam file from 10x Genomics RNA-Seq experiment and identify chromosome-wise pileups.

   The outputs of the pipeline contain a single 'Pyup.output.txt' of identified extreme pileup regions and quantification details for each chromosome.

   Edit config.yaml using a text editor to indicate input file location, chromosome of interests, and to change parameters.

   Note: the chromesomal names should be the same as in the header of input bam, i.e. "chr1" is considered different from "1".


Dependencies
------------
The pipeline is dependent on `conda3`, `snakemake`, `python3`, `bamtools`, `samtools`, `umi_tools`
