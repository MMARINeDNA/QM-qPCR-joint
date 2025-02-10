# Generate the Bash script for primer trimming
scripttext <- paste0(
  "#!/bin/bash\n\n",
  "PRIMERSEQ_F=", PRIMERSEQ_F, "\n",
  "PRIMERSEQ_R=", PRIMERSEQ_R, "\n",
  "PRIMERNAME=", PRIMERNAME, "\n\n",
  "CUTADAPT=", shQuote(CUTADAPT), "\n",
  "FASTQ_LOCATION=", shQuote(FASTQ_LOCATION), "\n\n",

  "cd \"${FASTQ_LOCATION}\"\n",
  "FILES=$(ls ./${PRIMERNAME}'-'*R1*.fastq.gz)\n",
  "mkdir '../for_dada2/'\n",
  "mkdir '../logs/'\n\n",

  "for i in $FILES\n",
  "do\n",
  "    FILE_NAME=$(echo ${i} | cut -d _ -f 1,2,3) # Grab everything before 'R1' (sample name)\n",
  "    echo $FILE_NAME\n",
  "    R1=${i}\n",
  "    R2=$(echo ${i} | sed 's/R1/R2/g')\n",
  "    ${CUTADAPT} -g ${PRIMERSEQ_F} \\\n",
  "               -G ${PRIMERSEQ_R} \\\n",
  "               -o ../for_dada2/${R1} \\\n",
  "               -p ../for_dada2/${R2} \\\n",
  "               --discard-untrimmed \\\n",
  "               -j 0 \\\n",
  "               ${R1} ${R2} 1>> '../logs/cutadapt_trim_report.txt'\n",
  "done\n"
)

# Write out the script to a file in code_etc
cat(
  scripttext,
  file = paste0(CODE_LOCATION, "/trim_primers.sh")
)
