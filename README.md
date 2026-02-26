# Comp-Bio-Pipeline-ProjectFASTQ files were retrieve via prefetch and fasterq-dump. SRR's were retrieved using prefetch {insert SRR} -O . in order to download SRR's into the current folder. SRA files were converted to FASTQ files with the command fasterq-dump --split-3 {Insert SRR}/{Insert SRR.sra}
SRRs were retrieved using the command prefetch {insert SRR} -O . in order to download SRRs to current folder. Then FASTQ files were obtained via fasterq-dumo --split-3 {insert SRR/insert SRR/sra}.
Dependencies: Kallisto, Sleuth, bowtie2, spades, biopython, blastn

To run the pipeline on the test samples

git clone https://github.com/jcapecci09/Comp-Bio-Pipeline-Project.git
cd Comp-Bio-Pipeline-Project
snakemake --cores 4 --config sample_folder=samples
