import argparse
import subprocess
import gzip
import os
import pandas as pd
import more_itertools

def clean_data(input_file, bed_file, chunksize=3_000_000):
    """
        Modifies Position and Locus from Neale Lab summary statistics using liftover-mapped positions, processing in chunks.
        Removes unnecessary BED file columns before loading to save memory.
    """

    print(f"Lifting over file: {input_file} using BED mapping from: {bed_file}")
    base_name = os.path.basename(input_file).replace(".vcf.gz", "")
    output_file = f"{base_name}_lifted.vcf.gz"
    bf = pd.read_csv(
        bed_file, sep='\t', header=None,
        usecols=[2, 3, 4],  # Load only "pos38", "chr", "old_pos"
        names=['pos38', 'chr', 'old_pos'],
        dtype={'pos38': 'int32', 'chr': 'category', 'old_pos': 'int32'}
    )
    bf['locus'] = bf['chr'].astype(str) + ':' + bf['old_pos'].astype(str)
    locus_to_pos38 = dict(zip(bf['locus'], bf['pos38']))
    del bf  # Free memory

    with gzip.open(input_file, mode='rt') as fp:
        fp = more_itertools.peekable(fp)
        header_line = next(fp).strip().split('\t')
    chunk_iterator = pd.read_csv(input_file, sep='\t', compression='gzip', chunksize=chunksize)

    for chunk_num, df in enumerate(chunk_iterator):
        print(f"Processing chunk {chunk_num + 1}...")

        # Assign correct column names
        df.columns = header_line

        # Create 'locus' column
        df.rename(columns={'pos': 'pos37'}, inplace=True)

        # Map 'pos37' to new 'pos' using 'locus' as the key
        df['pos'] = df['locus'].map(locus_to_pos38).fillna(999999999).astype("Int64")
        # Drop the old 'locus' column
        df.drop(columns=['locus'], inplace=True)
        # Create a new 'locus' column with 'chr' and the new 'pos'
        df['locus'] = df['chr'].astype(str) + ':' + df['pos'].astype(str)

        # Write the cleaned chunk to disk
        write_mode = 'wt' if chunk_num == 0 else 'at'
        write_header = chunk_num == 0  # Header only for the first chunk

        with gzip.open(output_file, write_mode) as f_out:
            df.to_csv(f_out, sep='\t', index=False, header=write_header)

        del df  # Free memory

    print(f"Final cleaned file saved as: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Clean Neale Lab Summary Statistics using Liftover")
    parser.add_argument('input_file', type=str, help="Path to the Neale Lab summary statistics file")
    parser.add_argument('bed_file', type=str, help="Path to the liftover BED file")

    args = parser.parse_args()
    clean_data(args.input_file, args.bed_file)


if __name__ == "__main__":
    main()