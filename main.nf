nextflow.enable.dsl=2

params.host_fasta = "${baseDir}/ref/hg38.fa"
params.kraken_db = "/Users/srinath/Documents/Folders/Python/Python_projects/meta_nextflow/db"
params.accessions = "SRR37617232"
params.output_dir = "results"

sra_ch = Channel.fromList(params.accessions.tokenize(","))

process DOWNLOAD_DATA {
    tag "Downloading ${sra_id}"
    
    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}*.fastq"), emit: raw_reads

    script:
    """
    # 1. Create a fake home environment inside the work directory
    export HOME=\$PWD
    mkdir -p \$HOME/.ncbi

    # 2. Write the config file (this is what fixed the Status 78 error!)
    echo '/LIBS/GUID = "8fa7e565-14f7-418f-9a91-450f3c647b74"' > \$HOME/.ncbi/user-settings.mkfg
    echo '/config/default = "true"' >> \$HOME/.ncbi/user-settings.mkfg

    # 3. Use the older fastq-dump because it supports the -X limit flag
    fastq-dump ${sra_id} --split-files -X 1000000
    """
}

process RUN_FASTQC{
    tag "FASTQ on ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(reads), emit: raw_reads
    path "*.html"

    script:
    """
    fastqc ${reads[0]} ${reads[1]} -o .
    """
}

process RUN_FASTP {
    tag "Trimming ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/trimmed_res", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    // These 'trimmed' reads are the ones we use for the rest of the pipeline
    tuple val(sample_id), path("trimmed_{1,2}.fastq"), emit: trimmed_reads
    path "*.{html,json}"

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
          -o trimmed_1.fastq -O trimmed_2.fastq \
          --html ${sample_id}_fastp_report.html \
          --json ${sample_id}_fastp_report.json
    """
}

process PREPARE_REFERENCE {
    tag "Building or Finding Index"
    
    output:
    // This 'emit' name must match exactly what you call in the workflow block
    path "hg38_index*", emit: bowtie2_index

    script:
    """
    # 1. Check if the FASTA file exists in your local project ref folder
    if [ -f "/Users/srinath/Documents/Folders/Python/Python_projects/meta_nextflow/ref/hg38.fa" ]; then
        echo "Found local FASTA, copying..."
        cp /Users/srinath/Documents/Folders/Python/Python_projects/meta_nextflow/ref/hg38.fa hg38.fa
    else
        echo "FASTA not found, downloading..."
        wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        mv Homo_sapiens.GRCh38.dna.primary_assembly.fa hg38.fa
    fi
    
    # 2. Build the index (Bowtie2 is smart; if the index exists, it will overwrite or you can add check logic)
    bowtie2-build --threads ${task.cpus} hg38.fa hg38_index
    """
}

process DEHOST_READS {
    tag "Dehosting ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/dehosted", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index_files // This receives the 'host_index' channel from above

    output:
    // We only care about the reads that DID NOT map (the microbes)
    tuple val(sample_id), path("dehosted_*.fastq.gz"), emit: dehosted_reads

    script:
    """
    # -x hg38_index: The prefix of your index files
    # --un-conc-gz: Saves the unmapped (microbial) paired reads to compressed fastq
    bowtie2 -x hg38_index \
            -1 ${reads[0]} -2 ${reads[1]} \
            --very-sensitive-local \
            --threads ${task.cpus} \
            --un-conc-gz dehosted_%.fastq.gz \
            -S /dev/null

    """
}

process TAXONOMIC_CLASSIFICATION {
    tag "Classifying ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/taxonomy", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path db // Passes the database directory into the process

    output:
    path "${sample_id}_kraken_report.txt", emit: report
    // path "${sample_id}_krona.html", emit: html

    script:
    """
    # 1. Run Kraken2
    # We use --use-names to get scientific names instead of just TaxIDs
    kraken2 --db ${db} \
            --paired ${reads[0]} ${reads[1]} \
            --threads ${task.cpus} \
            --use-names \
            --report ${sample_id}_kraken_report.txt \
            --output ${sample_id}_kraken_output.out

    # 2. Prepare the Krona input
    # This extracts the percentage and the name for the pie chart
    # awk -F'\t' '{print \$1 "\t" \$6}' ${sample_id}_kraken_report.txt > krona_input.txt

    # 3. Generate the Interactive Krona Chart
    # ktImportText krona_input.txt -o ${sample_id}_krona.html
    """
}

process ASSEMBLE_METAGENOME {
    tag "Assembling ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("metaspades_assembly/scaffolds.fasta"), emit: scaffolds
    path "metaspades_assembly/*"

    script:
    """
    # -m 8: Limits memory to 8GB (adjust based on your Mac's RAM)
    # -k 21,33,55: SPAdes uses different 'k-mer' lengths to build the assembly
    metaspades.py -1 ${reads[0]} -2 ${reads[1]} \
                  -o metaspades_assembly \
                  -t ${task.cpus} \
                  -m 8
    """
}

process EVALUATE_ASSEMBLY {
    tag "MetaQUAST on ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/quast_report", mode: 'copy'

    input:
    tuple val(sample_id), path(scaffolds)

    output:
    path "quast_results/*"

    script:
    """
    # --max-ref-number 10: Tells MetaQUAST to find the top 10 closest reference genomes
    # to compare your assembly against.
    metaquast.py ${scaffolds} \
                -o quast_results \
                --threads ${task.cpus} \
                --max-ref-number 10
    """
}

process FUNCTIONAL_ANNOTATION {
    tag "Annotating ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/annotation", mode: 'copy'

    input:
    tuple val(sample_id), path(scaffolds)

    output:
    path "functional_res/*"

    script:
    """
    prokka ${scaffolds} \
           --outdir functional_res \
           --prefix ${sample_id}_anno \
           --metagenome \
           --cpus ${task.cpus} \
           --compliant \
           --force \
           --usegenus
    """
}

process ALIGN_READS {
    tag "Aligning ${sample_id}"
    
    input:
    tuple val(sample_id), path(scaffolds)
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("mapped.sam"), emit: sam_file

    script:
    """
    bowtie2-build --threads ${task.cpus} ${scaffolds} scaffolds_index
    bowtie2 -x scaffolds_index -1 ${reads[0]} -2 ${reads[1]} --threads ${task.cpus} -S mapped.sam
    """
}

process CONVERT_AND_SORT_BAM {
    tag "Processing BAM for ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/mapping", mode: 'copy'

    input:
    tuple val(sample_id), path(sam)

    output:
    // We add the .stats file here so MultiQC can read it
    tuple val(sample_id), path("mapped_sorted.bam"), path("mapped_sorted.bam.bai"), emit: bam_files
    path "*.stats", emit: stats

    script:
    """
    samtools view -bS ${sam} | samtools sort -o mapped_sorted.bam
    samtools index mapped_sorted.bam
    samtools flagstat mapped_sorted.bam > ${sample_id}_mapping.stats
    """
}

process BINNING_METABAT2 {
    tag "Binning ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/bins", mode: 'copy'

    input:
    tuple val(sample_id), path(scaffolds)
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("bins_out/*.fa"), emit: bins

    script:
    """
    # 1. Summarize depth (how many reads are on each scaffold)
    jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam}

    # 2. Run MetaBAT2
    mkdir bins_out
    metabat2 -i ${scaffolds} -a depth.txt -o bins_out/bin
    """
}

process CHECKM_QUALITY {
    tag "Checking ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/checkm", mode: 'copy'

    input:
    tuple val(sample_id), path(bins)

    output:
    path "checkm_results/*"

    script:
    """
    mkdir bin_dir
    cp ${bins} bin_dir/
    
    # Use taxonomy_wf instead of lineage_wf
    # 'domain Bacteria' is much lighter on RAM
    checkm taxonomy_wf domain Bacteria bin_dir/ checkm_results/ -t ${task.cpus} -x fa
    """
}


process RUN_MULTIQC {
    publishDir "${params.output_dir}/multiqc", mode: 'copy'

    input:
    // .collect() gathers reports from all samples into one list
    path qc_files 

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

workflow {
    // 1. Create the starting channel from your SRA list
    sra_ch = Channel.fromList(params.accessions.tokenize(','))

    // 2. Start the pipeline
    DOWNLOAD_DATA(sra_ch)
    
   // Connect DOWNLOAD_DATA output to FASTQC
    // We use .out[0] because the tuple is the first thing emitted
    RUN_FASTQC(DOWNLOAD_DATA.out[0])
    
    // Connect FASTQC output to FASTP
    RUN_FASTP(RUN_FASTQC.out[0])
    
    // 4. Host Removal (Logic: Needs reads + reference index)
    PREPARE_REFERENCE()
    DEHOST_READS( 
        RUN_FASTP.out.trimmed_reads, 
        PREPARE_REFERENCE.out.bowtie2_index.collect() 
    )
    
    // 5. Taxonomy (Who is there?)
    TAXONOMIC_CLASSIFICATION(DEHOST_READS.out.dehosted_reads, params.kraken_db)
    
    // 6. Assembly & Evaluation
    ASSEMBLE_METAGENOME(DEHOST_READS.out.dehosted_reads)
    EVALUATE_ASSEMBLY(ASSEMBLE_METAGENOME.out.scaffolds)
    
    // 7. Functional Annotation (What are they doing?)
    FUNCTIONAL_ANNOTATION(ASSEMBLE_METAGENOME.out.scaffolds)
    
    // 8. Binning Prep (Split into two steps)
    ALIGN_READS(
        ASSEMBLE_METAGENOME.out.scaffolds, 
        DEHOST_READS.out.dehosted_reads
    )
    
    CONVERT_AND_SORT_BAM(
        ALIGN_READS.out.sam_file
    )

    // Update the input for Binning to use the new BAM files
    BINNING_METABAT2(
        ASSEMBLE_METAGENOME.out.scaffolds, 
        CONVERT_AND_SORT_BAM.out.bam_files
    )

    // 9. Final Quality Check on Bins
    CHECKM_QUALITY(BINNING_METABAT2.out.bins)

    // 10. MultiQC
    all_reports = RUN_FASTQC.out[1]
    .mix(RUN_FASTP.out[1])
    .mix(EVALUATE_ASSEMBLY.out)
    .mix(TAXONOMIC_CLASSIFICATION.out.report)
    // Use .map{ it[1] } to grab just the BAM/BAI files and skip the Sample ID
    .mix(CONVERT_AND_SORT_BAM.out.bam_files.map{ it[1] }) 
    .mix(CONVERT_AND_SORT_BAM.out.stats)
    .flatten()
    .collect()

    RUN_MULTIQC(all_reports)
}