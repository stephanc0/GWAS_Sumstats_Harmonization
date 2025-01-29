import argparse
import subprocess
import gzip
import os
import pandas as pd


def clean_data(input_file, bed_file):
    # Load the input file and bed file
    df = pd.read_csv(input_file, sep='\t', compression='gzip')
    bf = pd.read_csv(bed_file, sep='\t', header=None, names=['chromosome', 'new_posm1', 'pos38', 'chr', 'old_pos'])

    # Create 'locus' column in  bf
    bf['locus'] = bf['chr'].astype(str) + ':' + bf['old_pos'].astype(str)

    # Map the 'pos38' from `bf` to `df` based on matching 'locus'
    df['pos38'] = df['locus'].map(bf.set_index('locus')['pos38'])
    df['pos38'] = df['pos38'].fillna(-1).astype(int)
    df['pos38'] = df['pos38'].replace(-1, pd.NA)

    base_name = os.path.basename(input_file).split('.vcf.gz')[0]
    final_output_file = f"{base_name}_lifted.vcf.gz"
    df.to_csv(final_output_file, sep='\t', index=False, compression='gzip')

    print(f"Cleaned file saved as: {final_output_file}")

    return df, final_output_file

def main():
    parser = argparse.ArgumentParser(description="Clean Neale Lab Summary Statistics")
    parser.add_argument('input_file', type=str, help="Path to the Neale Lab summary statistics file")
    parser.add_argument('bed_file', type=str, help="Path to the Neale Lab summary statistics file")

    args = parser.parse_args()

    clean_data(args.input_file, args.bed_file)

if __name__ == "__main__":
    main()