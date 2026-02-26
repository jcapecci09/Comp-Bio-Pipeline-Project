# Define samples, their conditions and the genome accession code for this analysis
samples = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
conditions = {'SRR5660030':2, 'SRR5660033':6, 'SRR5660044':2, "SRR5660045":6}
genome_accession = "GCF_000845245.1"
IN_FOLDER = config.get("samples_folder", "samples") 

# Rule all to run pipeline all at once
# PipelineReport.txt is the expected final output
rule all:
    input:
        'PipelineReport.txt'

# Download the assembly from NCBI
rule download_assembly:
    output:
        "assembly.zip"
    shell:
        "datasets download genome accession {genome_accession} --include gff3,rna,cds,protein,genome,seq-report --filename {output}"

# unzip the assembly and put it in the assembly folder
rule unzip_assembly:
    input: 
        'assembly.zip'
    output:
        directory("assembly")
    shell:
        "unzip -o {input} -d {output}"

# Count the number of cds in assembly and create new fasta files with all transcripts
rule counting_cds:
    input:
        "assembly"
    output:
        tran = 'transcriptome.fasta',
        cds = 'counting_cds.txt'
    shell:
        "python Counting_cds.py -i {input}/ncbi_dataset/data/{genome_accession}/cds_from_genomic.fna -o {output.tran} -f {output.cds}"

# Build reference index for kallisto
rule kallisto_index:
    input:
        'transcriptome.fasta'
    output:
        'output.idx'
    shell:
        'kallisto index -i {output} {input}'

# Quantify number of transcripts per million
rule quantify_tpm:
    input:
        index = 'output.idx',
        fastq1 = IN_FOLDER + '/{sample}_1.fastq',
        fastq2 = IN_FOLDER + '/{sample}_2.fastq'
    output:
        directory('results/{sample}')
    shell:
        'time kallisto quant -i {input.index} -o {output} -b 30 -t 4 {input.fastq1} {input.fastq2}'

# Create table that will be used by sleuth to test for differential expression
rule create_table:
    input:
        sample_paths = expand('results/{sample}', sample=samples)
    output:
        'table.txt'
    shell:
        """
        echo -e 'sample\tcondition\tpath' > {output}
        echo -e 'SRR5660030\t2\tresults/SRR5660030' >> {output}
        echo -e 'SRR5660033\t6\tresults/SRR5660033' >> {output}
        echo -e 'SRR5660044\t2\tresults/SRR5660044' >> {output}
        echo -e 'SRR5660045\t6\tresults/SRR5660045' >> {output}
        """

# Calculate differential expression
rule perform_sleuth: 
    input:
        'table.txt'
    output:
        'sleuth_results.txt'
    shell:
        'Rscript differential_expression.R'

# Build index for bowtie
rule build_index:
    input:
        "assembly"
    output:
        "hcvm.1.bt2",
        "hcvm.2.bt2",
        "hcvm.3.bt2",
        "hcvm.4.bt2",
        "hcvm.rev.1.bt2",
        "hcvm.rev.2.bt2"
    shell:
        """
        bowtie2-build \
        assembly/ncbi_dataset/data/{genome_accession}/{genome_accession}_ViralProj14559_genomic.fna \
        hcvm
        """
# perform bowtie
rule map_reads:
    input:
        index="hcvm.1.bt2",
        read1= IN_FOLDER + "/{sample}_1.fastq",
        read2= IN_FOLDER + "/{sample}_2.fastq"
    output:
        sam="alignments/HCVMmap_{sample}.sam",
        bam="alignments/{sample}_mapped.bam",
        counts="alignments/read_counts_{sample}.txt"
    shell:
        """
        mkdir -p alignments
        before=$(($(wc -l < {input.read1}) / 4))

        bowtie2 --quiet -x hcvm -1 {input.read1} -2 {input.read2} -S {output.sam}

        samtools view -b -F 4 -f 2 {output.sam} > {output.bam}

        after=$(($(samtools view -c {output.bam}) / 2))

        echo "Sample {wildcards.sample}: $before pairs before, $after mapped reads." > {output.counts}
        """

# Aggregate read counts into one file
rule aggregate_counts:
    input:
        expand("alignments/read_counts_{sample}.txt", sample=samples)
    output:
        'read_counts.txt'
    shell:
        "cat {input} > {output}"

# BAM to FASTQ
rule bam_to_fastq:
    input:
        bam = "alignments/{sample}_mapped.bam"
    output:
        f1 = "fastq_for_spades/{sample}_R1.fastq",
        f2 = "fastq_for_spades/{sample}_R2.fastq"
    shell:
        """
        mkdir -p fastq_for_spades
        samtools sort -n {input.bam} | \
        samtools fastq -1 {output.f1} -2 {output.f2} -0 /dev/null -s /dev/null -n
        """
# Perform spades
rule spades:
    input:
        r1= IN_FOLDER + "/{sample}_1.fastq",
        r2= IN_FOLDER + "/{sample}_2.fastq"
    output:
        "spades_output/{sample}/contigs.fasta"
    shell:
        """
        spades.py \
            -1 {input.r1} \
            -2 {input.r2} \
            -o spades_output/{wildcards.sample} \
            -k 127 \
            -t 2 \
            --only-assembler
        """

# Find the longest contig from spades output
rule longest_contigs:
    input:
        contigs="spades_output/{x}/contigs.fasta"
    output:
        longest="spades_output/{x}/longest_contig.fasta"
    shell:
        "python largest_contig.py -i {input.contigs} -o {output.longest}"

# Obtain the database sequences in a fasta 
rule obtain_database_sequences:
    input:
        expand("spades_output/{x}/longest_contig.fasta", x=samples)
    output:
        "Betaherpesvirinae.fasta"
    shell:
        """
        python database.py -o {output}
        """
        
# Make database
rule make_database:
    input:
        'Betaherpesvirinae.fasta'
    
    output:
        "Betaherpesvirinae.nin",
        "Betaherpesvirinae.ndb",
        "Betaherpesvirinae.nhr",
        "Betaherpesvirinae.not",
        "Betaherpesvirinae.nsq",
        "Betaherpesvirinae.ntf",
        "Betaherpesvirinae.nto"
    
    shell:
        """
        makeblastdb -in {input} -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl
        """

# Run blast using the database 
rule run_blast:
    input:
        'Betaherpesvirinae.nin',
        fasta = "spades_output/{sample}/longest_contig.fasta"
    output:
        "blast_output/{sample}.tsv"
    shell:
        """
        blastn \
        -query {input.fasta} \
        -db Betaherpesvirinae \
        -out {output} \
        -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" \
        -max_hsps 1 \
        -max_target_seqs 5
        """
# aggrevate tsvs from blast output
rule aggregate_tsvs:
    input:
        csvs = expand("blast_output/{sample}.tsv", sample=samples)
    output:
        "all_blast_results.tsv"
    shell:
        "cat {input.csvs} > {output}"

# Add headers to blast output
rule add_headers:
    input:
        'all_blast_results.tsv'
    output:
        'all_blast_results.txt'
    shell:
        'python add_headers.py -i {input} -o {output}'

# Create final report
rule create_report:
    input:
        cds = 'counting_cds.txt',
        sleuth = 'sleuth_results.txt',
        counts = 'read_counts.txt',
        blast = 'all_blast_results.txt'
    output:
        'PipelineReport.txt'
    shell:
        "cat {input.cds} {input.sleuth} {input.counts} {input.blast} > {output}"
