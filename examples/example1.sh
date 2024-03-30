#!/bin/bash

output_dir="../output"
template_file="../template_example.fasta"

if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

python3 ../fastq_generator.py \
        --output_file "${output_dir}/example1.fastq" \
        --template_sequence $template_file \
        --with_parameters \
        --upstream_sequence "GCCGGCCATGGCG" \
        --left_margin 15 \
        --total_length 500 \
        --number_of_sequences 200 \
        --center 150 \
        --min_val 40 \
        --max_val 73 \
        --std_dev 75 \
        --noise_level 0.1 

        