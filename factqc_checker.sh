
#!/bin/bash
# Create output folder (safe if it already exists)
mkdir -p /home/ade/Renemma_2/CleanData/fastqc-output
set -Eeuo pipefail
# Run FastQC on the two inputs, extract PNGs/metrics folders
fastqc --extract -t 12 -o /home/ade/Renemma_2/CleanData/fastqc-output \
  /home/ade/Renemma_2/CleanData/LCS10337_PR1_Clean_Data1.fq.gz \
  /home/ade/Renemma_2/CleanData/LCS10337_PR1_Clean_Data2.fq.gz

# Resulting files (per sample) in fastqc-output/:
#   <sample>_fastqc.html
#   <sample>_fastqc.zip
#   <sample>_fastqc/Images/*.png   # because of --extract
#   <sample>_fastqc/fastqc_data.txt, summary.txt


# assuming you wrote FastQC outputs to this folder: following up with these will combine all the folders into one easy to ready quality check
multiqc /home/ade/Renemma_2/CleanData/fastqc-output \
  -o /home/ade/Renemma_2/CleanData/fastqc-output/summary
