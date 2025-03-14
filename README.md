This repository can be used to streamline and liftover summary statistics from various Biobanks.

Finngen summary statistics can be downloaded and processed like this.

Presumeably, this code will work on future Finngen data freezes if the column names are maintained.  View
column names in wsl terminal with 

wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_BMI_IRN.gz

zcat finngen_R12_BMI_IRN.gz | head -n 10

If downloaded using the gsutil link, there will be additional columns, which can be ignored.

IEUopenGWAS contains many GWAS from a variety of biobanks or consortiumns. 

Neale Lab and MRCIEU Summary Statistics can both be downloaded and processed like this:

wget https://gwas.mrcieu.ac.uk/files/ukb-a-61/ukb-a-61.vcf.gz

zcat ukb-a-61.vcf.gz | head -n 200

Run the clean_ukb.py like this: 

python clean_ukb.py ukb-a-61.vcf.gz --prevalence=0.26

Liftover can be run like this:

bash download.sh

bash run_liftover.sh ukb-a-61_cleaned.vcf.gz

Then to move over the new coordinates:

python liftover_ukb.py ukb-a-61_cleaned.vcf.gz ukb-a-61_cleaned_hg38.bed

Cleaned summary statistics files can be premunged for easy munging and analysis by LDSC:

python premunge_cleaned.py ukb-a-61_cleaned_lifted.vcf.gz --sample_size=335000





