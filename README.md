# Metagenomic Analysis Pipeline (Nextflow DSL2)

## Overview

This repository contains a **production-style metagenomic analysis pipeline** built using Nextflow (DSL2). The pipeline performs a complete workflow starting from raw sequencing data to taxonomic profiling, assembly, binning, and functional annotation.

It is designed to be:

* Modular
* Reproducible (Docker-based)
* Scalable (but currently configured for local execution)

---

## ⚠️ Important Note

Due to **local system storage and resource limitations**, this pipeline:

* Processes only **1,000,000 reads (1M reads)** per sample using `fastq-dump -X`
* Is configured for **sequential execution (queueSize = 1)**

This makes it suitable for **testing, development, and demonstration**, not full-scale production runs.

---

## Workflow Summary

### 1. Data Download

* Tool: `sra-tools`
* Downloads reads from SRA using accession IDs
* Limited to **1M reads**

### 2. Quality Control

* Tool: `FastQC`
* Generates quality reports for raw reads

### 3. Read Trimming

* Tool: `fastp`
* Removes low-quality bases and adapters

### 4. Host Removal (Dehosting)

* Tool: `Bowtie2`
* Aligns reads to human reference genome (hg38)
* Keeps **unmapped reads (microbial reads)**

### 5. Taxonomic Classification

* Tool: `Kraken2`
* Identifies microbial composition

### 6. Assembly

* Tool: `metaSPAdes`
* Assembles reads into contigs/scaffolds

### 7. Assembly Evaluation

* Tool: `MetaQUAST`
* Evaluates assembly quality

### 8. Functional Annotation

* Tool: `Prokka`
* Annotates genes and features

### 9. Read Mapping to Assembly

* Tool: `Bowtie2 + Samtools`
* Maps reads back to contigs

### 10. Binning

* Tool: `MetaBAT2`
* Groups contigs into genome bins

### 11. Bin Quality Assessment

* Tool: `CheckM`
* Evaluates completeness and contamination

### 12. Reporting

* Tool: `MultiQC`
* Aggregates all reports

---

## Reference & Database Setup

### Human Reference Genome (hg38)

The pipeline uses **hg38** for host removal.

#### Auto-download (handled in pipeline)

Downloads from Ensembl:

```
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

---

### Kraken2 Database

You must manually download and provide a Kraken2 database.

#### Recommended (Standard DB)

```bash
kraken2-build --standard --db kraken_db
```

#### Alternative (Prebuilt DB)

Download from:
https://benlangmead.github.io/aws-indexes/k2

#### Configure path in pipeline:

```groovy
params.kraken_db = "/path/to/kraken_db"
```

---

## Requirements

### Software

* Nextflow (>= 24.x)
* Docker (enabled)

### System Requirements

* macOS/Linux
* Minimum:

  * 8 GB RAM
  * 20–50 GB free storage (even with 1M reads)

---

## Running the Pipeline

```bash
nextflow run main.nf
```

---

## Configuration

### Key Parameters

```groovy
params.accessions = "SRR37617232"
params.output_dir = "results"
params.kraken_db = "/path/to/db"
params.host_fasta = "ref/hg38.fa"
```

---

## Output Structure

```
results/
 └── SRRXXXXX/
      ├── fastqc/
      ├── trimmed_res/
      ├── dehosted/
      ├── taxonomy/
      ├── assembly/
      ├── quast_report/
      ├── annotation/
      ├── mapping/
      ├── bins/
      └── checkm/
 └── multiqc/
```

---

## Execution Model

* Executor: `local`
* Queue: **1 job at a time**
* Docker containers used for all processes
* Platform forced to:

```
linux/amd64
```

---

## Limitations

* Uses **subset of reads (1M)** → not biologically exhaustive
* Local execution only (no cloud/HPC yet)
* Kraken2 classification depends heavily on database quality
* Assembly quality limited due to reduced read depth

---

## Future Improvements

* Add support for:

  * Full dataset processing
  * Cloud/HPC execution
  * Workflow parameterization via config profiles
  * Krona visualization (currently commented)
  * Automated Kraken DB download

---

## Author

Developed as a **bioinformatics workflow project** demonstrating:

* Workflow orchestration (Nextflow DSL2)
* Tool integration
* Reproducible research practices

---

## License

MIT License 

---

