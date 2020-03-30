Bootstrap: docker
From: continuumio/miniconda3:4.8.2

%labels
    AUTHOR yh154

%environment
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This sets global environment variables for anything run within the container
    export PATH="/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:"
    unset CONDA_DEFAULT_ENV
    export ANACONDA_HOME=/opt/conda

%post
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is going to be executed after the base container has been downloaded
    export PATH=/opt/conda/bin:$PATH
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install --yes snakemake
    conda install --yes bamtools=2.5.1
    conda install --yes samtools=1.9
    conda install --yes umi_tools=1.0.1
    conda install --yes pymc3=3.8
    conda clean --index-cache --tarballs --packages --yes

%runscript
    exec "$@"
