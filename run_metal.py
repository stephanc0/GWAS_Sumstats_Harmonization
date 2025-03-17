import argparse
import os
import subprocess
import gzip
import pandas as pd
import more_itertools


def preprocess_file(input_file, marker, chunksize=3_000_000):
    """
    Reads a GWAS summary statistics file in chunks and processes it according to METAL requirements.
    - If `marker == rsid`: removes rows with missing `rsid`, and drops `pos` and `locus`.
    - If `marker == locus`: removes rows with missing `locus`, converts `locus` to string, and drops `rsid` and `pos`.
    """
    base_name, ext = os.path.splitext(input_file)  # Splits ".gz" from "ukb-b-10215_cleaned_lifted.vcf.gz"
    base_name, main_ext = os.path.splitext(base_name)  # Splits ".vcf" from "ukb-b-10215_cleaned_lifted"
    output_file = f"{base_name}_metal{main_ext}{ext}"

    with gzip.open(input_file, mode='rt') as fp:
        fp = more_itertools.peekable(fp)
        header_line = next(fp).strip().split('\t')  # Read first line as header

    chunk_iterator = pd.read_csv(input_file, sep='\t', compression='gzip', chunksize=chunksize)

    for chunk_num, df in enumerate(chunk_iterator):
        print(f"Processing chunk {chunk_num + 1} for {input_file}...")

        df.columns = header_line  # Assign column names

        if marker == "rsid":
            # Drop rows where `rsid` is missing
            df.drop(
                df[
                    df['rsid'].isna() |
                    (df['rsid'] == ".") |
                    (df['rsid'] == "") |
                    (df['rsid'] == "rs999999999")
                    ].index,
                inplace=True
            )
            df.drop(columns=['pos', 'locus'], inplace=True, errors='ignore')

        elif marker == "locus":
            # Drop rows where `locus` is missing
            df = df[df['locus'].notna() & ~df['locus'].str.endswith(":999999999")]
            df['locus'] = df['locus'].astype(str)  # Ensure it's a string
            df.drop(columns=['rsid', 'pos'], inplace=True, errors='ignore')

        # Write processed chunk to a temporary file
        write_mode = 'wt' if chunk_num == 0 else 'at'
        write_header = chunk_num == 0  # Write header only for the first chunk

        with gzip.open(output_file, write_mode) as f_out:
            df.to_csv(f_out, sep='\t', index=False, header=write_header)

        del df  # Free memory

    print(f"Preprocessed file saved as: {output_file}")
    return output_file  # Return the name of the modified file


def generate_metal_script(file1, file2, marker, output_file, metal_exec):
    """
    Generates a METAL script to perform meta-analysis on two GWAS summary statistics files.
    """
    metal_script = f"""
    SCHEME STDERR
    MARKER {marker}
    ALLELE ref alt
    EFFECT beta
    STDERR sebeta
    PVAL pval
    COLUMNCOUNTING LENIENT
    PROCESS {file1}
    PROCESS {file2}
    OUTFILE {output_file} .tbl
    ANALYZE
    QUIT
    """
    script_filename = "metal_script.txt"
    with open(script_filename, "w") as f:
        f.write(metal_script)

    print(f"METAL script written to {script_filename}")

    # Run METAL
    print("Running METAL...")
    subprocess.run(f"{metal_exec} < {script_filename}", shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Run METAL meta-analysis for two GWAS summary statistics files.")
    parser.add_argument("file1", type=str, help="First GWAS summary statistics file")
    parser.add_argument("file2", type=str, help="Second GWAS summary statistics file")
    parser.add_argument("marker", type=str, choices=["rsid", "locus"],
                        help="Column to use as the marker for matching variants")
    parser.add_argument("output_file", type=str, help="Base name for the output file")
    parser.add_argument("--metal_exec", type=str, default="./metal/metal",
                        help="Path to METAL executable (default: ./metal/metal)")

    args = parser.parse_args()

    # Preprocess both input files
    cleaned_file1 = preprocess_file(args.file1, args.marker)
    cleaned_file2 = preprocess_file(args.file2, args.marker)

    # Run METAL on the cleaned files
    generate_metal_script(cleaned_file1, cleaned_file2, args.marker, args.output_file, args.metal_exec)

    # Remove temporary files
    os.remove(cleaned_file1)
    os.remove(cleaned_file2)
    print("Temporary files deleted.")


if __name__ == "__main__":
    main()

