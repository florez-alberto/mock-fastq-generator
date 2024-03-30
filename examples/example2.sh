#!/bin/bash

output_dir="../output"
template_file="../template_example.fasta"
parameters_file="../parameters_example.json"

if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

python3 ../fastq_generator.py \
        --output_file "${output_dir}/example2.fastq" \
        --template_sequence $template_file \
        --parameters_file $parameters_file

        