#!/usr/bin/env python3

import csv
import os
import re

# file splitting functions
def split_fastq_path_into_name_components(fastq_file):
    bn = os.path.basename(fastq_file)
    name_components = bn.rsplit(".")[0].split("_")
    return name_components


# globals
read_dir = 'data/reads'
sample_key = 'data/sample_key.csv'

# find fastq files
read_dir_files = [(dirpath, filenames) for (dirpath, dirnames, filenames)
                  in os.walk(read_dir)]
fastq_files = []
for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            fastq_files.append(os.path.join(dirpath, filename))

# get a list of all samples
name_components = {x: split_fastq_path_into_name_components(x)
                   for x in fastq_files}
all_samples = set(name_components[x][0] for x in name_components)

# get the set of barcodes per sample
sample_to_bc = {}
for sample in all_samples:
    sample_to_bc[sample] = set(name_components[x][1]
        for x in name_components
        if name_components[x][0] == sample)

# get the set of lanes per sample
sample_to_lane = {}
for sample in all_samples:
    sample_to_lane[sample] = set(name_components[x][2]
        for x in name_components
        if name_components[x][0] == sample)



# generate a dict of sample-to-name
with open(sample_key) as csvfile:
    csvreader = csv.reader(csvfile, delimiter = ",")
    next(csvreader)
    sample_to_name = {x: y for x, y in csvreader}


rule all:
    input: 
        expand('output/merged/{all_samples}_R1.fastq',
            all_samples = all_samples),
        expand('output/merged/{all_samples}_R2.fastq',
            all_samples = all_samples)

# I think I will need a dicts of sample : bc, sample : lane, sample : fn

rule merge_per_sample:
    input:
        r1 = 'data/reads/{{all_samples}}_{bc}_{lane}_R1_{fn}.fastq.gz',
        r2 = 'data/reads/{{all_samples}}_{bc}_{lane}_R2_{fn}.fastq.gz'
    output:
        r1 = touch('output/merged/{all_samples}_R1.fastq'),
        r2 = touch('output/merged/{all_samples}_R2.fastq')
    shell:
        'echo "'
        'zcat {input.r1} >'
        '{output.r1} ; '
        'zcat {input.r2} >'
        '{output.r2}"'

