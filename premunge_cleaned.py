import argparse
import gzip
import os
import pandas as pd
from scipy.stats import norm
import more_itertools

def process_chunk(df, sample_size):
    """
    Processes a chunk of the GWAS summary statistics file.
    - Computes Z-scores
    - Computes p-values if missing
    - Filters valid SNPs
    - Renames columns to match LDSC format
    """
    # Calculate Z scores
    df['Z'] = (df['beta'] / df['sebeta']).round(7)

    # Calculate p-values if they are not present
    if 'pval' not in df.columns:
        df['pval'] = 2 * norm.sf(abs(df['Z']))
        df['pval'] = df['pval'].apply(lambda x: '{:.7e}'.format(x) if x < 0.0001 else round(x, 7))

    # Ensure p-values are within valid range [0, 1]
    df = df[(df['pval'].astype(float) >= 0) & (df['pval'].astype(float) <= 1)]

    # Add sample size
    df['N'] = sample_size

    # Define valid bases for A1 and A2
    valid_bases = {'A', 'C', 'G', 'T'}

    # Filter valid SNPs
    df = df[(df['alt'].isin(valid_bases)) & (df['ref'].isin(valid_bases))]

    # Rename columns to match LDSC format
    df = df[['rsid', 'alt', 'ref', 'af_alt', 'beta', 'pval', 'N', 'Z']].rename(
        columns={'rsid': 'SNP', 'alt': 'A1', 'ref': 'A2', 'af_alt': 'MAF', 'beta': 'BETA', 'pval': 'P'}
    )
    # Filter out variants with invalid RSIDS
    df = df[df['SNP'] != "rs999999999"]

    return df


def premunge_data(input_file, sample_size, chunksize=3_000_000):
    """
    Reads the input VCF file in chunks, processes each chunk, and writes to a final munged file.
    This avoids splitting the file into separate smaller files.
    """
    base_name = os.path.basename(input_file).replace(".vcf.gz", "")
    output_file = f"{base_name}_premunged.vcf.gz"

    # Use peekable to extract the header first
    with gzip.open(input_file, mode='rt') as fp:
        fp = more_itertools.peekable(fp)
        header_line = next(fp).strip().split('\t')  # Extract header

    chunk_iterator = pd.read_csv(input_file, sep='\t', compression='gzip', chunksize=chunksize)

    first_chunk = True
    with gzip.open(output_file, 'wt') as f_out:
        for chunk_num, df in enumerate(chunk_iterator):
            print(f"Processing chunk {chunk_num+1}...")

            df.columns = header_line  # Assign proper column names

            df = process_chunk(df, sample_size)  # Process the data

            # Write the processed chunk to disk
            df.to_csv(f_out, sep='\t', index=False, header=first_chunk)

            first_chunk = False  # Ensure only the first chunk writes the header

            del df  # Free memory

    print(f"Final premunged file saved as: {output_file}")
    return output_file


def main():
    parser = argparse.ArgumentParser(description="Premunge large GWAS summary statistics files in chunks")
    parser.add_argument('input_file', type=str, help="Path to input VCF file")
    parser.add_argument('--sample_size', type=int, required=True, help="Sample size for the GWAS")

    args = parser.parse_args()

    # Munge data using chunk processing
    premunge_data(args.input_file, args.sample_size)


if __name__ == "__main__":
    main()
