This directory contains a metabarcoding pipeline that runs from demultiplexed FASTQ files through primer trimming, DADA2 analysis, and local BLAST-based taxonomic assignment. 
Below is an overview of each file and its function.

1. make_primer_shell_script.R
Purpose
	Generates a bash script (trim_primers.sh) to trim primer sequences using Cutadapt.
How It Works
	Reads user-defined primer sequences (forward & reverse), the location of cutadapt, and the FASTQ folder path.
	Produces a shell script that loops through matching FASTQ files, trimming primers and discarding untrimmed reads.
Usage
	Called automatically by the main pipeline script (Pipeline_code.R) to ensure reproducible primer trimming.

2. Pipeline_code.R
Purpose
	The master script of the pipeline that orchestrates:
	Primer trimming (via trim_primers.sh).
	Quality filtering & ASV inference (DADA2).
	Hash key generation and ASV table creation.
	Local BLAST of novel sequences and taxonomic assignment (via LCA_function.R).
	Producing final taxon/haplotype tables.
Key Steps
	Input Preparation: Loads primer data (primer.data.csv), sets user-specific paths, and copies necessary scripts to a temporary code directory.
	DADA2 Processing: Filters, denoises, merges, and removes chimeras from sequences.
	Local BLAST (via Local_blast_function.R): Searches sequences against a local BLAST database.
	LCA Analysis (via LCA_function.R): Assigns sequences to taxa using a lowest common ancestor approach.
	Output: Generates final taxonomic tables in both long and wide formats, plus optional haplotype-level tables.
Usage
	Ensure you have Blast, taxonkit, and cutadapt available globally in your machine. I recomment running their respective tutorials first
	Update paths to your local environment (e.g., path to FASTQ files, cutadapt, taxids_to_exclude.txt if used, etc.).
	Run in R (≥ 4.0) with packages: tidyverse, dada2, digest, seqinr, ShortRead, here, etc.

3. primer.data.csv
Purpose
	Holds primer metadata used to trim and filter sequences.
Contents
	name: Unique primer name (e.g., MFU).
	seq_f, seq_r: Forward & reverse primer sequences.
	primer_length_f, primer_length_r: Length of each primer.
	max_amplicon_length, min_amplicon_length, overlap: Amplicon constraints for DADA2 merging.
Usage
	Pipeline_code.R and make_primer_shell_script.R reference these fields to trim and filter reads consistently.
	Ensure the PRIMERNAME you set in Pipeline_code.R matches the name column here.

4. taxids_to_exclude.txt
Purpose
	An optional file listing taxonomic IDs to exclude during BLAST.
Contents
	One taxid per line (commented lines starting with # or blank lines are ignored).
Usage
	Read by Local_blast_function.R (or the older ceg_blast() approach) to pass -negative_taxids to BLAST.
	If not needed, set this file to NULL in Pipeline_code.R.

5. LCA_function.R
Purpose
	Implements the Lowest Common Ancestor (LCA) algorithm to assign sequences to taxa after BLAST.
Key Steps
	Reads the BLAST output (in outfmt 6 format).
	Uses a taxonomic database (e.g., NCBI) to find the lowest shared rank for each query.
	Updates the master database (_database.csv) with newly assigned sequences.
Usage
	Called by Pipeline_code.R after the local BLAST step to finalize taxonomic designations for each ASV.

6. Local_blast_function.R
	Purpose
	Provides a local BLAST function (local_blast()) that replaces remote SSH-based BLAST calls.
	Functionality
	Reads negative taxids from taxids_to_exclude.txt (if provided) to create a comma-separated list.
	Invokes blastn on your local machine using user-defined parameters (e.g., perc_identity, word_size, evalue, etc.).
	Writes a standard BLAST output file (e.g., new_annotations.txt) in outfmt 6.
	Usage
	Ensure you have BLAST+ installed locally and that the local BLAST database (e.g., nt_euk) is available.
	Called within Pipeline_code.R to annotate newly discovered ASVs.

Additional Notes
Dependencies
	R packages: dada2, tidyverse, digest, seqinr, ShortRead, here, etc.
	External tools: cutadapt, BLAST+, taxonkit. Make sure you have them runnig on the global environment
Running the Pipeline
	Modify Pipeline_code.R to reflect your absolute or relative paths.
	Run Pipeline_code.R in an R session (≥ 4.0).
	Inspect logs and outputs in the generated logs and outputs folders.