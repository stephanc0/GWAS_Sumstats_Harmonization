This repository can be used to streamline, liftover, and meta_analyze summary statistics from Finngen, as well as the UK Biobank, 
curated either by Neale Lab, or by MRC-IEU. These summary statistics can be downloaded, respectively, with commands 
like this:

wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_BMI_IRN.gz

wget https://gwas.mrcieu.ac.uk/files/ukb-a-61/ukb-a-61.vcf.gz

Summary Statistics can be streamlined(unnecessary columns removed and important statistics in their own columns) like: 

python clean_finngen.py finngen_R12_BMI_IRN.gz --sd=.14
python clean_ukb.py ukb-a-61.vcf.gz --prevalence=0.26

UKBB files use the GRCh37 coordinates, and can be lifted over. First, download the necessary files.  Next, liftover the
coordinates using a .bed file as required.  Then, convert the coordinates in your summary statistic file.

bash download.sh

bash run_liftover.sh ukb-a-61_cleaned.vcf.gz

python liftover_ukb.py ukb-a-61_cleaned.vcf.gz ukb-a-61_cleaned_hg38.bed

Cleaned summary statistics files can be premunged for easy munging and analysis by LDSC:

python premunge_cleaned.py ukb-a-61_cleaned_lifted.vcf.gz --sample_size=335000

Cleaned summary statistics files can also be meta-analyzed with eachother, using either locus (chr:pos) or using 
rsid:

python run_metal.py ukb-b-10215_cleaned_lifted.vcf.gz finngen_R12_BMI_IRN_cleaned.tsv.gz rsid HGSWTMETA







