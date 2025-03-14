import argparse
import subprocess
import gzip
import os
import pandas as pd

def clean_data(input_file, p):
    temp_file = "temp_cleaned.vcf.gz"
    command = f"zcat {input_file} | grep -v '^##' | gzip > {temp_file}"
    subprocess.run(command, shell=True, check=True)
    os.replace(temp_file, input_file)

    df = pd.read_csv(input_file, sep='\t')
    assert len(df.FORMAT.unique()) == 1

    df_split = df.iloc[:, -1].str.split(':', expand=True)
    df_split.columns = ['beta', 'sebeta', 'mlogp', 'af_alt', 'sample_size', 'id']
    df = pd.concat([df, df_split[['beta', 'sebeta', 'mlogp', 'af_alt', 'sample_size']]], axis=1)
    # print("Unique IDs in 'sample sizes' column:", df['sample_size'].unique())
    # print("Number of unique sample sizes:", df['sample_size'].nunique())
    assert df['sample_size'].nunique() == 1

    df['beta'] = pd.to_numeric(df['beta'])
    df['sebeta'] = pd.to_numeric(df['sebeta'])
    df['mlogp'] = pd.to_numeric(df['mlogp'])
    df['af_alt'] = pd.to_numeric(df['af_alt'])

    df.drop(columns=['QUAL', 'INFO', 'FILTER', 'FORMAT', 'UKB-a-61', 'sample_size'], inplace=True)
    df.rename(columns={'#CHROM': 'chr', 'POS': 'pos', 'REF': 'ref', 'ALT': 'alt',}, inplace=True)
    df['locus'] = df['chr'].astype(str) + ':' + df['pos'].astype(str)

    lin_to_log_conversion = p * (1-p)

    df.rename(columns={'beta': 'lin_beta'}, inplace=True)
    df.rename(columns={'sebeta': 'lin_sebeta'}, inplace=True)
    df['beta'] = df['lin_beta'] / lin_to_log_conversion
    df['sebeta'] = df['lin_sebeta'] / lin_to_log_conversion
    df.drop(columns=['lin_beta', 'lin_sebeta'], inplace=True)

    base_name = os.path.basename(input_file).split('.vcf.gz')[0]
    final_output_file = f"{base_name}_cleaned.vcf.gz"
    df.to_csv(final_output_file, sep='\t', index=False, compression='gzip')

    print(f"Cleaned file saved as: {final_output_file}")

    return df, final_output_file

def main():
    parser = argparse.ArgumentParser(description="Clean Neale Lab Summary Statistics")
    parser.add_argument('input_file', type=str, help="Path to the Neale Lab summary statistics file")
    parser.add_argument('--prevalence', type=float, required=True, help="Prevalence of the trait in the population")

    args = parser.parse_args()

    clean_data(args.input_file, args.prevalence)

if __name__ == "__main__":
    main()