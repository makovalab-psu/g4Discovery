import argparse
import os
import pandas as pd

from g4DiscoveryFuncs import *
from g4HunterFuncs import *

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fa", "--fasta_file", type=str, required=True, 
                        help="Path to the input FASTA file")
    parser.add_argument("-chr", "--chromosome", type=str, required=True, 
                        help="Chromosome identifier, either an integer or a single letter")
    parser.add_argument("-o", "--output", type=str, required=True, 
                        help="Path to the output BED file")
    parser.add_argument("-t", "--tetrad", type=int, default=3, required=False,
                        help="Minimum number of tetrads for a G4 to be considered")
    parser.add_argument("-ps", "--pqsscore", type=int, default=40, required=False,
                        help="Minimum pqsfinder score for a G4 to be considered")
    parser.add_argument("-hs", "--g4hunter", type=float, default=1.5, required=False,
                        help="Minimum absolute G4Hunter score for a G4 to be considered")
    parser.add_argument("-psd", "--docker_min_pqsscore", type=str, required=False, default=30,
                        help="Minimum pqsfinder score for the docker to run")
    args = parser.parse_args()

    input_file_path = args.fasta_file
    output_file_path = os.path.dirname(args.output)
    output_file_name = os.path.basename(args.fasta_file) + ".pqs"

    print("Running Docker container: kxk302/pqsfinder:1.0.0")
    run_docker(input_file_path, output_file_path, output_file_name, pqs_min_score=args.docker_min_pqsscore)

    print(f'Using the following parameters: tetrad={args.tetrad}, pqsscore={args.pqsscore}, g4hunter={args.g4hunter}')

    # Validate chr argument to ensure it is either an integer or a single letter
    if not (args.chromosome.isdigit() or (len(args.chromosome) == 1 and args.chromosome.isalpha())):
        parser.error("The -chr argument must be either an integer or a single letter.")

    # Filter G4s
    filteredG4s = filterG4s(fasta_file=os.path.join(output_file_path, output_file_name), chr=args.chromosome, min_tetrad=args.tetrad, min_score=args.pqsscore, min_g4hunterscore=args.g4hunter)

    # Filter the DataFrames
    plus_strand_df = filterNonOverlappingG4sCmplx(filteredG4s[filteredG4s["strand"] == "+"])
    minus_strand_df = filterNonOverlappingG4sCmplx(filteredG4s[filteredG4s["strand"] == "-"])

    # Write the output to a file if there are G4s found on either strand
    if plus_strand_df.empty and minus_strand_df.empty:
        print("No G4s found on either strand, based on the filtering criteria.")
    else:
        final = pd.concat([df for df in [plus_strand_df, minus_strand_df] if not df.empty], axis=0)
        # Sort the final dataframe by start and end positions
        final.sort_values(by=["start","end"], ascending=[True,False], inplace=True)
        final.to_csv(args.output, sep="\t", header=False, index=False) 
