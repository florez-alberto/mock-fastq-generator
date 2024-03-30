
import argparse
import utils as utils

# import difflib
# template_sequence= "TGGATAACGCACCCTATTGTGGGTCGATGTGGAGTACGAATCGAAAAAATTAACCAGAGCCTG"
# seqs=utils.add_noise_to_template_sequence(template_sequence, noise_level=0.1)
# print(template_sequence)
# print(seqs)
# # Assuming template_sequence and seqs are your two strings
# differences = difflib.ndiff(template_sequence, seqs)
# print(''.join(differences))
# exit()
def main(args):
    import utils as utils
    import json
    
    if args.template_sequence:
        template_sequence= utils.read_fasta_file(args.template_sequence)
    else:
        template_sequence= utils.read_fasta_file("./template_example.fasta")

    if args.with_parameters:
        upstream_sequence = args.upstream_sequence
        left_margin= args.left_margin
        total_length= args.total_length
        number_of_sequences= args.number_of_sequences
        center = args.center
        min_val = args.min_val
        max_val = args.max_val
        std_dev = args.std_dev
        noise_level = args.noise_level
    else:
        with open(args.parameters_file) as f:
            parameters = json.load(f)
            upstream_sequence = parameters['upstream_sequence']
            left_margin= parameters['left_margin']
            total_length= parameters['total_length']
            number_of_sequences= parameters['number_of_sequences']
            center = parameters['center']
            min_val = parameters['min_val']
            max_val = parameters['max_val']
            std_dev = parameters['std_dev']
            noise_level = parameters['noise_level']
    
    assembled_sequences = utils.assemble_sequences(template_sequence,left_margin,upstream_sequence,total_length, number_of_sequences, center, min_val, max_val, std_dev,noise_level) 
    utils.write_fastq_file(args.output_file, assembled_sequences)
    return

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="\n Generate a fastq file with sequences assembled from a template sequence and an upstream sequence.")
    argparser.add_argument("--output_file", type=str, default="output.fastq", help="Output file")
    argparser.add_argument("--template_sequence", type=str, default=None, help="Template sequence in fasta format")

    argparser.add_argument("--with_parameters", action='store_true', help="Use the parameters in the arguments")
    argparser.add_argument("--upstream_sequence", type=str, default="GCCGGCCATGGCG")
    argparser.add_argument("--left_margin", type=int, default=15)
    argparser.add_argument("--total_length", type=int, default=500)
    argparser.add_argument("--number_of_sequences", type=int, default=10)
    argparser.add_argument("--center", type=int, default=150)
    argparser.add_argument("--min_val", type=int, default=40)
    argparser.add_argument("--max_val", type=int, default=73)
    argparser.add_argument("--std_dev", type=int, default=75)
    argparser.add_argument("--noise_level", type=float, default=0.1)
    argparser.add_argument("--parameters_file", type=str, default="parameters.json", help="Parameters file")

    args = argparser.parse_args()

    if not args.with_parameters:
            print("The --with_parameters flag was not passed. The conditional arguments will not be used.")

    main(args)  
