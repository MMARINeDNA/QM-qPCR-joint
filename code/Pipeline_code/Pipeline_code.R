##############################################################################
##                         Welcome to The Kelly Lab's                       ##
##              Metabarcoding Taxonomic Assignment Pipeline                 ##
## This pipeline will take you from de-multiplexed fastq files to a matrix  ##
## with taxon names and counts. It should work for Linux or macOS.          ##
## PC users may look into WSL2 to run a similar environment.                ##
##############################################################################
##    Before you start, ensure you have access to the following files:      ##
## 1. Pipeline_code.R - This very code (you are here, so you got it!)       ##
## 2. Your Fastqs - Names must start with the primer (e.g. MFU_123.fastq.gz)##
## 3. make_primer_shell_script.R - A script that writes another script      ##
## 4. primer_data.csv - Primer info (sequences, lengths, etc.)              ##
## 5. Local_blast_function.R - Custom function to run BLAST locally         ##
## 7. LCA_function.R - Function to perform Lowest Common Ancestor analysis  ##
## 8. taxids_to_exclude.txt (optional) - List of taxids to exclude    	    ##
##                                                                          ##
## Required Programs:                                   					##
## 1. Cutadapt: https://cutadapt.readthedocs.io/en/stable/installation.html ## 
## 2. Blast:    https://www.ncbi.nlm.nih.gov/books/NBK569861/               ##
## 3. Taxonkit: https://bioinf.shenwei.me/taxonkit/                         ##
##############################################################################

# Restart your R session before running each time to avoid environment carryover.
# .rs.restartR()  # If running in RStudio
# Be cautious when clearing your environment:
# rm(list = ls())

# Mark this file so 'here' package knows our relative location.
here::i_am("pipeline_code.R")

# Load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(dada2))
suppressMessages(library(digest))
suppressMessages(library(seqinr))
suppressMessages(library(ssh))
suppressMessages(library(sys))
suppressMessages(library(ShortRead))
suppressMessages(library(here))

################################################################################
# 1. Define user-specific paths and parameters
################################################################################

# Point to your code location, databases, and other paths:
# (Replace "/path/to/..." with actual paths on your system.)

# The directory containing a subfolder called "Fastq" with your .fastq or .fastq.gz files
PARENT_LOCATION <- "/path/to/parent_directory/"  # e.g. "/mnt/c/Users/YourName/Pipeline/"
if (strsplit(PARENT_LOCATION, "")[[1]][nchar(PARENT_LOCATION)] != "/") {
  PARENT_LOCATION <- paste0(PARENT_LOCATION, "/")
}

# A unique name for this run
RUN_NAME <- "Example_run_name"
# The primer name must match an entry in your primer_data.csv file
PRIMERNAME <- "MFU"

# Threshold sequence similarity for annotation (e.g. 96 for species-level match)
SP_THOLD <- 96

# Optional file listing taxids to exclude from BLAST
NEGATIVE_TAXIDS_FILE <- "taxids_to_exclude.txt"  # or set to NULL if not used

# Path to a local or shared database for storing/reusing annotated sequences
DATABASE_LOCATION <- "/path/to/databases/"  # e.g. "/mnt/c/Users/YourName/Pipeline/databases/"
if (strsplit(DATABASE_LOCATION, "")[[1]][nchar(DATABASE_LOCATION)] != "/") {
  DATABASE_LOCATION <- paste0(DATABASE_LOCATION, "/")
}

# Final directory for processed outputs
PROCESSED_LOCATION <- "/path/to/processed_output/"  # e.g. "/mnt/c/Users/YourName/Pipeline/processed/"

# Create an output folder for this run
system2("mkdir", shQuote(paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME)))

# Subfolders used by the pipeline
FASTQ_LOCATION    <- paste0(PARENT_LOCATION, "Fastq")     # Contains your raw fastq files
CODE_LOCATION     <- paste0(PARENT_LOCATION, "code_etc")  # Folder for code and scripts
TRIMMED_LOCATION  <- paste0(PARENT_LOCATION, "for_dada2") # Folder after primer trimming
FILTERED_LOCATION <- paste0(PARENT_LOCATION, "filtered")   # Folder for quality-filtered files
OUTPUT_LOCATION   <- paste0(PARENT_LOCATION, "outputs")    # Folder for pipeline outputs

# Create needed subfolders
system2("mkdir", args=shQuote(CODE_LOCATION))
system2("mkdir", args=shQuote(TRIMMED_LOCATION))
system2("mkdir", args=shQuote(FILTERED_LOCATION))
system2("mkdir", args=shQuote(OUTPUT_LOCATION))

# Path to cutadapt and taxonkit executables
CUTADAPT <- "/usr/local/bin/cutadapt"  # Adjust if installed elsewhere
TAXONKIT_PATH <- "/usr/local/bin/taxonkit"  # Adjust if installed elsewhere

# Check if cutadapt and taxonkit are installed
system2(CUTADAPT, args="--version")
system2(TAXONKIT_PATH, args="version")

#calling custom  for blast and LCA
source(here("Local_blast_function.R"))
source("LCA_function.R")

################################################################################
# 2. Load primer data and set primer-specific parameters
################################################################################

# Load primer data (ensure primer_data.csv is accessible)
primer.data <- read.csv("/path/to/primer_data.csv")

# Extract relevant info for the chosen primer
PRIMERSEQ_F <- primer.data %>% filter(name == PRIMERNAME) %>% pull(seq_f)
PRIMERSEQ_R <- primer.data %>% filter(name == PRIMERNAME) %>% pull(seq_r)
PRIMERLENGTH_F <- primer.data %>% filter(name == PRIMERNAME) %>% pull(primer_length_f)
PRIMERLENGTH_R <- primer.data %>% filter(name == PRIMERNAME) %>% pull(primer_length_r)
MAX_AMPLICON_LENGTH <- primer.data %>% filter(name == PRIMERNAME) %>% pull(max_amplicon_length)
MIN_AMPLICON_LENGTH <- primer.data %>% filter(name == PRIMERNAME) %>% pull(min_amplicon_length)
OVERLAP <- primer.data %>% filter(name == PRIMERNAME) %>% pull(overlap)

# Copy key scripts and primer data to code_etc folder
write.csv(primer.data, paste0(CODE_LOCATION, "/primer_data.csv"))
system2("cp", args=c(shQuote(here("make_primer_shell_script.R")), shQuote(CODE_LOCATION)))
system2("cp", args=c(shQuote(here("pipeline_code_3.0.R")), shQuote(CODE_LOCATION)))
if (!is.null(NEGATIVE_TAXIDS_FILE)) {
  if (file.exists(here(NEGATIVE_TAXIDS_FILE))) {
    system2("cp", args=c(shQuote(here(NEGATIVE_TAXIDS_FILE)), shQuote(CODE_LOCATION)))
  }
}

# Source the script to generate the primer-trimming shell script
system2("cp", args=c(shQuote(here("make_primer_shell_script.R")), shQuote(CODE_LOCATION)))
source(paste0(CODE_LOCATION, "/make_primer_shell_script.R"))

# (Optional) If re-running on the same data, you can skip most of the pipelien and go straight to step 6:
SKIP_PIPE <- FALSE
if (SKIP_PIPE) {
  existing_asv_table <- file.path(PROCESSED_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME), "outputs", paste0(RUN_NAME, "_", PRIMERNAME, "_ASV_table.csv"))
  existing_hash_key  <- file.path(PROCESSED_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME), "outputs", paste0(RUN_NAME, "_", PRIMERNAME, "_hash_key.csv"))
  current_asv        <- read_csv(existing_asv_table)
  conv_table         <- read_csv(existing_hash_key)
}

################################################################################
# 3. Primer trimming with shell script
################################################################################

# The script "trim_primers.sh" is created by make_primer_shell_script.R
system2("sh", args=shQuote(paste0(CODE_LOCATION, "/trim_primers.sh")))

################################################################################
# 4. DADA2 pipeline: quality filtering, error learning, merging
################################################################################

# Get list of trimmed files
filelist <- system2("ls", args=shQuote(TRIMMED_LOCATION), stdout=TRUE)
fnFs <- filelist[str_detect(filelist, "_R1")]  # forward reads
fnRs <- filelist[str_detect(filelist, "_R2")]  # reverse reads

# Extract sample names
sample.names1 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names  <- paste(sample.names1, sample.names2, sep="_")

# Construct full paths for filtered files
filtFs <- file.path(FILTERED_LOCATION, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(FILTERED_LOCATION, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter out empty files (sometimes happens)
setwd(TRIMMED_LOCATION)
file.empty <- function(filenames) file.info(filenames)$size <= 20
empty_files <- file.empty(fnFs) | file.empty(fnRs)
fnFs <- fnFs[!empty_files]
fnRs <- fnRs[!empty_files]
filtFs <- filtFs[!empty_files]
filtRs <- filtRs[!empty_files]
sample.names <- sample.names[!empty_files]

# (Optional) Inspect quality profiles to decide truncation/trim
# plotQualityProfile(fnFs[1:3])
# plotQualityProfile(fnRs[1:3])

# filterAndTrim
out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen      = MIN_AMPLICON_LENGTH,
  maxN          = 0,
  maxEE         = c(2,2),
  truncQ        = 2,
  rm.phix       = TRUE,
  compress      = TRUE,
  multithread   = FALSE,
  matchIDs      = TRUE
)

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, selfConsist=TRUE, multithread=TRUE, MAX_CONSIST=20)
dadaRs <- dada(filtRs, err=errR, selfConsist=TRUE, multithread=TRUE, MAX_CONSIST=20)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap=OVERLAP, verbose=TRUE, trimOverhang=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
cat("Proportion of non-chimeric reads:", sum(seqtab.nochim) / sum(seqtab), "\n")

# Filter by size (optional)
indexes.to.keep <- which(
  nchar(colnames(seqtab.nochim)) <= MAX_AMPLICON_LENGTH &
  nchar(colnames(seqtab.nochim)) >= MIN_AMPLICON_LENGTH
)
cleaned.seqtab.nochim    <- seqtab.nochim[, indexes.to.keep]
filteredout.seqtab.nochim <- seqtab.nochim[, !indexes.to.keep]

# Save filtered-out table for reference
dir.create(paste0(PARENT_LOCATION, "logs"), showWarnings=FALSE)
write.csv(filteredout.seqtab.nochim, paste0(PARENT_LOCATION, "logs/filtered_out_asv.csv"))

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, paste0(PARENT_LOCATION, "logs/tracking_reads.csv"))

################################################################################
# 5. Hash key and ASV table creation
################################################################################
conv_file         <- file.path(OUTPUT_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME, "_hash_key.csv"))
conv_file.fasta   <- file.path(OUTPUT_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME, "_hash_key.fasta"))
ASV_file          <- file.path(OUTPUT_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME, "_ASV_table.csv"))

seqtab.nochim.df <- as.data.frame(cleaned.seqtab.nochim)

# Create a tibble with hash IDs for each unique sequence
conv_table <- tibble(
  Hash     = map_chr(colnames(seqtab.nochim.df), ~ digest(.x, algo="sha1", serialize=FALSE)),
  Sequence = colnames(seqtab.nochim.df)
)

# Write hash key
readr::write_csv(conv_table, conv_file)

# Write hash key to a FASTA
seqinr::write.fasta(
  sequences = as.list(conv_table$Sequence),
  names     = as.list(conv_table$Hash),
  file.out  = conv_file.fasta
)

# Create an ASV abundance table in long format
sample.df <- tibble::rownames_to_column(seqtab.nochim.df, "Sample_name")
sample.df <- data.frame(append(sample.df, c(Label=PRIMERNAME), after=1))

current_asv <- sample.df %>%
  pivot_longer(
    cols      = !c("Sample_name", "Label"),
    names_to  = "Sequence",
    values_to = "nReads"
  ) %>%
  filter(nReads > 0)

# Merge with the hash table
current_asv <- merge(current_asv, conv_table, by="Sequence") %>%
  dplyr::select(-Sequence) %>%
  dplyr::relocate(Hash, .after="Label")

# Write ASV table
readr::write_csv(current_asv, ASV_file)

################################################################################
# 6. BLAST-based taxonomic annotation (Local BLAST)
################################################################################

# Path to your local BLAST database (similar to BLAST_DB in your template)
LOCAL_BLAST_DB_PATH <- "/path/to/blastdb/nt_euk"  # Adjust as needed

# Check if there's a primer-specific database and make the list of new ASVs to be blasted
db_file <- paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv")

if (file.exists(db_file)) {
  db <- read.csv(db_file, row.names=1)
  seen    <- which(conv_table$Hash %in% db$Hash)
  notseen <- which(!conv_table$Hash %in% db$Hash)
  
  if (length(notseen) > 0) {
    seqinr::write.fasta(
      sequences = as.list(conv_table$Sequence[notseen]),
      names     = as.list(conv_table$Hash[notseen]),
      file.out  = paste0(OUTPUT_LOCATION, "/seqs_to_annotate.fasta")
    )
  } else {
    # If everything is already in the DB
    seqinr::write.fasta(
      sequences = as.list(conv_table$Sequence),
      names     = as.list(conv_table$Hash),
      file.out  = paste0(OUTPUT_LOCATION, "/seqs_to_annotate.fasta")
    )
  }
} else {
  # If no database exists, write all sequences
  seqinr::write.fasta(
    sequences = as.list(conv_table$Sequence),
    names     = as.list(conv_table$Hash),
    file.out  = paste0(OUTPUT_LOCATION, "/seqs_to_annotate.fasta")
  )
}


# Perform local BLAST on the new sequences
local_blast(
  PATH_TO_FASTA         = file.path(OUTPUT_LOCATION, "seqs_to_annotate.fasta"),
  DB_PATH               = LOCAL_BLAST_DB_PATH,
  OUTPUT_FILE           = file.path(OUTPUT_LOCATION, "new_annotations.txt"),
  NEGATIVE_TAXIDS_FILE  = if (!is.null(NEGATIVE_TAXIDS_FILE)) {
                            file.path(CODE_LOCATION, NEGATIVE_TAXIDS_FILE)
                          } else {
                            NULL
                          },
  PERCENT_IDENTITY      = 90,       # or 96 if you want to match e.g. vertebrate settings
  WORD_SIZE             = 30,
  EVALUE                = 1e-30,
  MAXIMUM_MATCHES       = 50,
  CULLING               = 50
)

# Now run LCA on the local BLAST results
annotation_file <- file.path(OUTPUT_LOCATION, "new_annotations.txt")
if (file.exists(annotation_file) && file.size(annotation_file) > 0) {
  if (file.exists(db_file)) {
    LCA(
      BLASTOUTPUT  = annotation_file,
      FASTA        = file.path(OUTPUT_LOCATION, "seqs_to_annotate.fasta"),
      DB_PATH_IN   = db_file,
      DB_PATH_OUT  = db_file,
      SP_THOLD     = SP_THOLD
    )
    db <- read.csv(db_file, row.names=1)
  } else {
    db <- LCA(
      BLASTOUTPUT = annotation_file,
      FASTA       = file.path(OUTPUT_LOCATION, "seqs_to_annotate.fasta"),
      SP_THOLD    = SP_THOLD
    )
    write.csv(db %>% distinct(), db_file)
  }
} else {
  if (file.exists(db_file)) {
    db <- read.csv(db_file, row.names=1)
  } else {
    stop("No database file found and annotation_file is empty.")
  }
}

################################################################################
# 7. Writing final taxon/haplotype tables
################################################################################

tax_table_file <- file.path(OUTPUT_LOCATION, paste0(RUN_NAME, "taxon_table.csv"))
tax_table_wide <- file.path(OUTPUT_LOCATION, paste0(RUN_NAME, "taxon_table_wide.csv"))

# Merge current_asv with database to get BestTaxon, Class, etc.
tax_table <- current_asv %>%
  left_join(db %>% dplyr::select(Hash, BestTaxon, Class), by="Hash") %>%
  group_by(Sample_name, BestTaxon, Class) %>%
  summarise(nReads = sum(nReads), .groups="drop")

# Write in long format
readr::write_csv(tax_table, tax_table_file)

# Write wide format
tax_table %>%
  tidyr::pivot_wider(
    id_cols     = c(BestTaxon, Class),
    names_from  = Sample_name,
    values_from = nReads,
    values_fill = 0
  ) %>%
  readr::write_csv(tax_table_wide)

# Also produce a haplotype-level table (if your db has HaplotypeNumber)
merged_data <- current_asv %>% left_join(db, by="Hash")
merged_data <- merged_data %>% mutate(BestTaxon_Haplotype = paste0(BestTaxon, "_", HaplotypeNumber))

haplotype_data <- merged_data %>%
  group_by(BestTaxon_Haplotype, Sample_name, Class) %>%
  summarise(nReads = sum(nReads, na.rm=TRUE), .groups="drop")

haplotype_table_path <- file.path(OUTPUT_LOCATION, paste0(RUN_NAME, "haplotype_table.csv"))
haplotype_data %>%
  tidyr::pivot_wider(
    id_cols     = c(BestTaxon_Haplotype, Class),
    names_from  = Sample_name,
    values_from = nReads,
    values_fill = 0
  ) %>%
  write.csv(haplotype_table_path, row.names=FALSE)

################################################################################
# 8. Cleanup (optional)
################################################################################
CLEANUP <- TRUE  # Toggle this off if you want to keep intermediate files

if (CLEANUP){
  for (folder in c("code_etc", "outputs", "filtered", "for_dada2", "logs")) {
    target <- paste0(PARENT_LOCATION, folder)
    if (dir.exists(target)) {
      system2("rm", args=c("-r", shQuote(target)))
    }
  }
}

message("Pipeline complete. Outputs (taxon/haplotype tables, hash keys, etc.) are in:\n", file.path(PROCESSED_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME)))
