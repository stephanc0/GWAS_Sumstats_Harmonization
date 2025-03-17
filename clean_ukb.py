import argparse
import more_itertools
import gzip
import os
import pandas as pd

def clean_data(input_file, p=None, sd=None, comment = b'##', chunksize=3_000_000):
    """
    Cleans the Neale Lab summary statistics and applies appropriate transformations in chunks.
    """
    base_name = os.path.basename(input_file).replace(".vcf.gz", "")
    output_file = f"{base_name}_cleaned.vcf.gz"
    capitalized_base_name = base_name[:3].upper() + base_name[3:]  # Capitalize first 3 letters

    with gzip.open(input_file, mode='r') as fp:
        fp = more_itertools.peekable(fp)

        while fp.peek().startswith(comment):
            next(fp)

        header_line = next(fp).decode().strip().split('\t') # .split('\t')
        chunk_iterator = pd.read_csv(input_file, sep='\t', chunksize=chunksize)

    for chunk_num, df in enumerate(chunk_iterator):
        print(f"Processing chunk {chunk_num+1}...")

        df.columns = header_line
        if 'FORMAT' not in df.columns:
            raise ValueError("FORMAT column missing from chunk, check header assignment")

        assert len(df.FORMAT.unique()) == 1

        format_fields = df['FORMAT'].iloc[0].split(":")        # Extracts ['ES', 'SE', 'LP', 'AF', 'ID']
        contains_SS = "SS" in format_fields  # Check if SS is in the FORMAT

        df_split = df[capitalized_base_name].str.split(':', expand=True)
        df = df.drop(columns=[capitalized_base_name])  # Remove original column

        if contains_SS:
            df_split.columns = ['beta', 'sebeta', 'mlogp', 'af_alt', 'sample_size', 'id']
            df = pd.concat([df, df_split[['beta', 'sebeta', 'mlogp', 'af_alt', 'sample_size']]], axis=1)
        else:
            df_split.columns = ['beta', 'sebeta', 'mlogp', 'af_alt', 'id']
            df = pd.concat([df, df_split[['beta', 'sebeta', 'mlogp', 'af_alt']]], axis=1)

        # Convert necessary columns to numeric
        numeric_cols = ['beta', 'sebeta', 'mlogp', 'af_alt']
        df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')

        # Drop unnecessary columns
        drop_cols = ['QUAL', 'INFO', 'FILTER', 'FORMAT']
        if contains_SS:
            drop_cols.append('sample_size')
        df.drop(columns=[col for col in drop_cols if col in df.columns], inplace=True)

        df.rename(columns={'#CHROM': 'chr', 'POS': 'pos', 'REF': 'ref', 'ALT': 'alt', 'ID':'rsid'}, inplace=True)
        df['locus'] = df['chr'].astype(str) + ':' + df['pos'].astype(str)

        # Apply logistic conversion if prevalence (p) is provided
        if p is not None:
            df.rename(columns={'beta': 'lin_beta', 'sebeta': 'lin_sebeta'}, inplace=True)
            lin_to_log_conversion = p * (1 - p)
            df['beta'] = df['lin_beta'] / lin_to_log_conversion
            df['sebeta'] = df['lin_sebeta'] / lin_to_log_conversion
            df.drop(columns=['lin_beta', 'lin_sebeta'], inplace=True)

        # Apply linear transformation if standard deviation (sd) is provided
        if sd is not None:
            df['beta'] *= sd
            df['sebeta'] *= sd

        # Write the cleaned chunk to disk
        write_mode = 'wt' if chunk_num == 0 else 'at'  # Overwrite first chunk, append others
        write_header = chunk_num == 0  # Write header only for the first chunk

        with gzip.open(output_file, write_mode) as f_out:
            df.to_csv(f_out, sep='\t', index=False, header=write_header)

        del df

    print(f"Final cleaned file saved as: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Clean Neale Lab Summary Statistics in chunks")
    parser.add_argument('input_file', type=str, help="Path to the Neale Lab summary statistics file")
    parser.add_argument('--prevalence', type=float, default=None, help="Prevalence of the trait (for logistic)")
    parser.add_argument('--sd', type=float, default=None, help="Standard deviation of the trait (for linear)")

    args = parser.parse_args()
    clean_data(args.input_file, p=args.prevalence, sd=args.sd)


if __name__ == "__main__":
    main()
