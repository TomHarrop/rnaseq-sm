#!/usr/bin/env python3

import csv
import os


#############
# FUNCTIONS #
#############

# file splitting functions
def split_fastq_path_into_name_components(fastq_file):
    bn = os.path.basename(fastq_file)
    name_components = bn.rsplit(".")[0].split("_")
    return name_components


# lookup wildcards for merge function
def merge_wildcard_resolver(wildcards):
    '''
    The merge_per_sample rule passes two wildcards, `sample_name` and `read`.
    Look up the file names using the name_to_sample dict and return a list of
    files that match by filename and read number
    '''
    my_sample = name_to_sample[wildcards.sample_name]
    matched_files = sorted(
        set(x for x in name_components
            if (name_components[x][0] == my_sample and
                name_components[x][3] == wildcards.read)))
    return matched_files

###########
# GLOBALS #
###########

read_dir = 'data/reads'
sample_key = 'data/sample_key.csv'
fasta_file = 'data/genome.fasta'
star_index = 'output/star_index/index'


########
# PREP #
########

# find fastq files
read_dir_files = [(dirpath, filenames) for (dirpath, dirnames, filenames)
                  in os.walk(read_dir)]
fastq_files = []
for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            fastq_files.append(os.path.join(dirpath, filename))

# split all filepaths into name components
name_components = {x: split_fastq_path_into_name_components(x)
                   for x in fastq_files}

# generate a dict of name-to-sample
# N.B. FAKE NAMES
with open(sample_key) as csvfile:
    csvreader = csv.reader(csvfile, delimiter = ",")
    next(csvreader)
    name_to_sample = {y: x for x, y in csvreader}

# get a list of all sample names
all_sample_names = list(set(x for x in name_to_sample))

#########
# RULES #
#########

rule all:
    input: 
        expand('output/star/{sample_name}.Aligned.sortedByCoord.out.bam',
               sample_name = all_sample_names)

rule index_genome:
    input:
        fasta_file
    output:
        touch(star_index)
    shell:
        'echo '
        'STAR '
        'genomeGenerate {input} '
        '> {output}'

rule merge_per_sample:
    input:
        merge_wildcard_resolver
    output:
        touch('output/merged/{sample_name}_{read}.fastq.gz')
    shell:
        'echo "'
        'cat {input} '
        '> {output}"'

rule bbduk:
    input:
        r1 = 'output/merged/{sample_name}_R1.fastq.gz',
        r2 = 'output/merged/{sample_name}_R2.fastq.gz'
    output:
        r1 = touch('output/bbduk/{sample_name}_R1.fastq.gz'),
        r2 = touch('output/bbduk/{sample_name}_R2.fastq.gz')
    shell:
        'echo '
        'bbduk '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        '"FIXME PARAMS"'

rule star:
    input:
        r1 = 'output/bbduk/{sample_name}_R1.fastq.gz',
        r2 = 'output/bbduk/{sample_name}_R2.fastq.gz',
        star_index = star_index
    output:
        touch('output/star/{sample_name}.Aligned.sortedByCoord.out.bam')
    shell:
        'echo '
        'STAR --genomeDir {input.star_index} '
        '--readFilesIn {input.r1} {input.r2} '
        '--outSAMtype BAM SortedByCoordinate '
        '--twopassMode Basic '
        '--runThreadN 6 '
        '--outFileNamePrefix {wildcards.sample_name}.'


