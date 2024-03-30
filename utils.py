import random
import numpy as np
import string

def read_fasta_file(path):
    NGS_sequences=[]

# open file and read the content in a list
    with open(path,"r") as filehandle:
        for line in filehandle:
            # remove linebreak which is the last character of the string
            currentPlace = line[:-1]
            # add item to the list
        NGS_sequences=currentPlace.upper()
    return NGS_sequences


def generate_phred_scores(num_points, center, min_val, max_val, std_dev,ASCII_BASE = 33):
    phred_scores_ascii = {i: chr(i + ASCII_BASE) for i in range(0, 74 - ASCII_BASE)}

    # Generate Gaussian curve
    x = np.linspace(0, num_points - 1, num_points)
    gaussian_curve = np.exp(-((x - center) ** 2) / (2 * std_dev ** 2))

    # Add noise
    noise = np.random.normal(0, 0.1, num_points)  # Adjust the second parameter for noise intensity
    smooth_noisy_curve = gaussian_curve + noise

    # Scale to desired range
    smooth_noisy_curve = smooth_noisy_curve * (max_val - min_val) + min_val

    # Ensure values are within range
    smooth_noisy_curve = np.clip(smooth_noisy_curve, min_val, max_val)

    #round to the nearest integer and convert to int
    smooth_noisy_curve = np.round(smooth_noisy_curve).astype(int)

    phred_scores = ''.join([phred_scores_ascii[value-ASCII_BASE] for value in smooth_noisy_curve ])

    return phred_scores

def generate_at_header():
    #make a random string @XXXX:1234
    return "@" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=4)) + ":"+ ''.join(random.choices(string.digits, k=4))

def add_noise_to_template_sequence(template_sequence, noise_level=0.1):
    bases= ['A', 'C', 'G', 'T']
    noisy_sequence = []
    for base in template_sequence:
        if random.random() < noise_level:
            noisy_sequence.append(random.choice(bases))
        else:
            noisy_sequence.append(base)
    return ''.join(noisy_sequence)

def assemble_sequences(template_sequence,left_margin=110,upstream_sequence="GCCGGCCATGGCG",total_length=500 ,number_of_sequences=10,center=150, min_val=40, max_val=73, std_dev=75, noise_level=0.1):
    bases= ['A', 'C', 'G', 'T']
    sequences = []
    for i in range(number_of_sequences):
        random_seq_left_margin = ''.join(random.choices(bases, k=left_margin))

        random_seq_right_margin = ''.join(random.choices(bases, k=total_length - left_margin - len(template_sequence+upstream_sequence)))

        noisy_template_sequence = add_noise_to_template_sequence(template_sequence, noise_level)

        sequence = upstream_sequence + random_seq_left_margin + noisy_template_sequence + random_seq_right_margin

        sequences.append(sequence)
    assembled_sequences=[]
    for sequence in sequences:
        assembler = []
        header = generate_at_header()
        assembler.append(header)
        assembler.append(sequence)
        assembler.append("+")
        total_length =len(sequence)
        phred_string=generate_phred_scores(total_length, center, min_val, max_val, std_dev)
        assembler.append(phred_string)
        assembled_sequences.append(assembler)

    return assembled_sequences

def write_fastq_file(file, sequences):
    with open(file, 'w') as f:
        for sequence in sequences:
            f.write(sequence[0]+'\n')
            f.write(sequence[1]+'\n')
            f.write(sequence[2]+'\n')
            f.write(sequence[3]+'\n')

