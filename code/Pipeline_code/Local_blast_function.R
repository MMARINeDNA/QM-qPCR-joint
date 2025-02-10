local_blast <- function(
    PATH_TO_FASTA,
    DB_PATH,
    OUTPUT_FILE = "new_annotations.txt",
    NEGATIVE_TAXIDS_FILE = NULL,
    PERCENT_IDENTITY = 90,
    WORD_SIZE = 30,
    EVALUE = 1e-30,
    MAXIMUM_MATCHES = 50,
    CULLING = 50
) {
  # Build negative_taxids option if file is provided
  negative_taxids_option <- ""
  if (!is.null(NEGATIVE_TAXIDS_FILE) && file.exists(NEGATIVE_TAXIDS_FILE)) {
    # Extract valid taxids (non-empty lines that are not commented)
    taxids <- readLines(NEGATIVE_TAXIDS_FILE, warn = FALSE)
    taxids <- gsub("#.*", "", taxids)               # Remove comments
    taxids <- taxids[nzchar(trimws(taxids))]        # Remove empty lines
    # Join them as a comma-separated list
    exclude_list <- paste(taxids, collapse = ",")
    # Build the BLAST argument
    negative_taxids_option <- paste0("-negative_taxids \"", exclude_list, "\"")
  }

  # Make sure EVALUE is a string for the system call
  EVALUE <- as.character(EVALUE)

  # Construct the system command for blastn
  blast_cmd <- paste(
    "blastn",
    "-query", shQuote(PATH_TO_FASTA),
    "-db", shQuote(DB_PATH),
    "-num_threads 16",
    "-perc_identity", shQuote(PERCENT_IDENTITY),
    "-word_size", shQuote(WORD_SIZE),
    "-evalue", shQuote(EVALUE),
    "-max_target_seqs", shQuote(MAXIMUM_MATCHES),
    "-culling_limit", shQuote(CULLING),
    negative_taxids_option,
    "-outfmt \"6 sscinames scomnames qseqid sseqid pident length mismatch gapopen",
    "qcovus qstart qend sstart send evalue bitscore staxids qlen qcovs\"",
    "-out", shQuote(OUTPUT_FILE)
  )

  message("Running local BLAST...\n", blast_cmd)
  # Run blastn via system2
  system(blast_cmd)
  message("Local BLAST complete. Results in: ", OUTPUT_FILE)
}
