import argparse
import gzip
import os
import pandas as pd
import more_itertools  # To handle header extraction consistently

def clean_data(input_file, chunksize=3_000_000):
    """
    Cleans FinnGen summary statistics file in chunks while ensuring correct header handling.
    """
    # Define output file name
    base_name = os.path.basename(input_file).replace(".gz", "")
    output_file = f"{base_name}_cleaned.tsv.gz"

    # Extract header line
    with gzip.open(input_file, mode='rt') as fp:
        fp = more_itertools.peekable(fp)
        header_line = next(fp).strip().split('\t')  # Read first line as header
    
    # Read the file in chunks **without skipping rows**
    chunk_iterator = pd.read_csv(input_file, sep='\t', compression='gzip', chunksize=chunksize)

    for chunk_num, df in enumerate(chunk_iterator):
        print(f"Processing chunk {chunk_num+1}...")

        # Assign column names manually (ensuring consistency across chunks)
        df.columns = header_line

        # Rename '#chrom' column to 'chr'
        if '#chrom' in df.columns:
            df.rename(columns={'#chrom': 'chr'}, inplace=True)

        # Drop 'nearest_genes' column if it exists
        if 'nearest_genes' in df.columns:
            df.drop(columns=['nearest_genes'], inplace=True)

        # Create a 'locus' column as 'chr:pos'
        df['locus'] = df['chr'].astype(str) + ':' + df['pos'].astype(str)
        df.rename(columns={'rsids': 'rsid'}, inplace=True)

        # Write chunk to output file (overwrite first chunk, append subsequent ones)
        write_mode = 'wt' if chunk_num == 0 else 'at'
        write_header = chunk_num == 0  # Write header only for the first chunk
        df['rsid'] = df['rsid'].fillna("rs999999999").replace(["", "."], "rs999999999")
        with gzip.open(output_file, write_mode) as f_out:
            df.to_csv(f_out, sep='\t', index=False, header=write_header)

        # Free memory
        del df

    print(f"Final cleaned file saved as: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Clean FinnGen Summary Statistics in chunks")
    parser.add_argument('input_file', type=str, help="Path to the FinnGen summary statistics file")

    args = parser.parse_args()
    clean_data(args.input_file)

if __name__ == "__main__":
    main()
