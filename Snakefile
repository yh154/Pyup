import os
import glob
from snakemake.utils import validate, min_version
from snakemake.shell import shell
from itertools import combinations


min_version("5.11.2")

configfile: "config.yaml"

CHROMS=config["chrom"].replace(" ", "")
SAMPLES=CHROMS.split(",")

rule all:
    input:
        "Pyup.output.txt"

rule samtools_split_bam:
    input:
        config["bam"]
    output:
        temp(touch("bamsplit.done"))
    conda:
        "envs/bamtools.yaml"
    threads:
        config["params"]["threads"]
    shell:
        "bamtools split -in {input} -stub \"chr\" -reference -refPrefix \"\""

rule samtools_filter:
    input:
        "bamsplit.done"
    output:
        header=temp("{sample}.header"),
        bam=temp("{sample}.filtered.bam")
    conda:
        "envs/samtools.yaml"
    params:
        "chr.{sample}.bam"
    shell:
        "samtools view -H {params} > {output.header} && "
        "samtools view -F 1284 -q 255 {params} | grep \"CB:Z:\" |"
        "grep \"UB:Z:\" | cat {output.header} - | samtools view -bS -o {output.bam} - && rm {params}"


rule samtools_sort:
    input:
        "{sample}.filtered.bam"
    output:
        temp("{sample}.sorted.bam")
    conda:
        "envs/samtools.yaml"
    params:
        "chr{sample}"
    shell:
        "samtools sort -T {wildcards.sample} -o {output} {input}"

rule samtools_index:
    input:
        "{sample}.sorted.bam"
    output:
        temp("{sample}.sorted.bam.bai")
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule umitools_dedup:
    input:
        idx="{sample}.sorted.bam.bai",
        bam="{sample}.sorted.bam"
    output:
        temp("{sample}.sorted.dedup.bam")
    conda:
        "envs/umitools.yaml"
    shell:
        "umi_tools dedup --extract-umi-method=tag --umi-tag=UB --cell-tag=CB -I {input.bam} -S {output}"

rule samtools_depth:
    input:
        "{sample}.sorted.dedup.bam"
    output:
        temp("{sample}.sorted.dedup.bam.nozero.depth")
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools depth {input} > {output}"

rule running_median:
    input:
        "{sample}.sorted.dedup.bam.nozero.depth"
    output:
        temp("chr.{sample}.roi.txt")
    conda:
        "envs/script.yaml"
    params:
        site_width=config["running_median"]["site_width"],
        ratio=config["running_median"]["ratio"],
        min_depth=config["running_median"]["min_depth"]
    script:
        "scripts/Pyup.py"

rule cat_roi:
    input:
        expand("chr.{sample}.roi.txt", sample=SAMPLES)
    output:
        "Pyup.output.txt"
    shell:
        "echo -e \"chr\tstart\tend\tid\tmedian_zscore\tstrand\twidth\tmedian_coverage\" | cat - {input} > {output}"
