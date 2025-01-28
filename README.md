This repository can be used to modify and streamline Summary statistics from various Biobanks.

Finngen summary statistics come with many metrics, some of which are useless or redundant. They can be viewed here: https://r12.finngen.fi/. 
They can be downloaded directly, with permission, like this example for BMI:

wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_BMI_IRN.gz

The columns are #chrom, pos, ref, alt, rsids, nearest_genes, pval, mlogp, beta, sebeta, and af_alt.  Presumeably, 
this code will work on future Finngen data freezes if the column names are maintained.  View
column names in wsl terminal with 

zcat finngen_R12_BMI_IRN.gz | head -n 10

If downloaded using the gsutil link, there will be additional columns, which can be ignored.

IEUopenGWAS contains many GWAS from a variety of biobanks or consortiumns.

For instance, Neale Lab GWAS can be downloaded like this:

wget https://gwas.mrcieu.ac.uk/files/ukb-a-61/ukb-a-61.vcf.gz

zcat ukb-a-61.vcf.gz | head -n 200

Run the script like this: 

python clean_neale_lab.py shortened_file.vcf.gz --prevalence 0.26