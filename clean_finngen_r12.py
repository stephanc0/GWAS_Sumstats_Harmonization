import argparse
import gzip
import pandas as pd
def clean_data(input_file):
    print(f"Cleaning the file: {input_file}")
    # Open and read the gzipped file
    with gzip.open(input_file, 'rt') as file:
        # Load the file into a pandas DataFrame
        # Assumes the file is tab-separated
        df = pd.read_csv(input_file, sep='\t')
    # Rename the '#chrom' column to 'chr'
    df.rename(columns={'#chrom': 'chr'}, inplace=True)

    # Remove the 'nearest_genes' column
    if 'nearest_genes' in df.columns:
        df.drop(columns=['nearest_genes'], inplace=True)

    # Create a 'locus' column as 'chr:pos'
    df['locus'] = df['chr'].astype(str) + ':' + df['pos'].astype(str)

    # Optionally, print the first few rows to verify the changes
    print(df.head())

    return df

def main():
    parser = argparse.ArgumentParser(description="Clean FinnGen Summary Statistics")
    parser.add_argument('input_file', type=str, help="Path to the FinnGen summary statistics file")

    args = parser.parse_args()
    clean_data(args.input_file)

if __name__ == "__main__":
    main()
