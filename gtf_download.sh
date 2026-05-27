#!/bin/bash

# Define the base directory
BASE_DIR="/home/ade/Renemma_2"

# Create the base directory if it does not exist
mkdir -p $BASE_DIR

# Download the GTF file for mRatBN7.2 from Ensembl release 113
wget -O $BASE_DIR/Rattus_norvegicus.mRatBN7.2.113.gtf.gz ftp://ftp.ensembl.org/pub/release-113/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.113.gtf.gz || { echo "Failed to download GTF file"; exit 1; }

# Decompress the GTF file
gunzip -c $BASE_DIR/Rattus_norvegicus.mRatBN7.2.113.gtf.gz > $BASE_DIR/Rattus_norvegicus.mRatBN7.2.113.gtf || { echo "Failed to decompress GTF file"; exit 1; }