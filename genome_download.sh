#!/bin/bash

# Define the base directory
BASE_DIR="/home/ade/Renemma_2"

# Create the base directory if it does not exist
mkdir -p $BASE_DIR

# Download rat genome files for mRatBN7.2 from Ensembl release 113
wget -O $BASE_DIR/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz ftp://ftp.ensembl.org/pub/release-113/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz || { echo "Failed to download cdna file"; exit 1; }
wget -O $BASE_DIR/Rattus_norvegicus.mRatBN7.2.ncrna.fa.gz ftp://ftp.ensembl.org/pub/release-113/fasta/rattus_norvegicus/ncrna/Rattus_norvegicus.mRatBN7.2.ncrna.fa.gz || { echo "Failed to download ncrna file"; exit 1; }
wget -O $BASE_DIR/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz ftp://ftp.ensembl.org/pub/release-113/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz || { echo "Failed to download dna file"; exit 1; }

# Create decoy list for Salmon
grep "^>" <(gunzip -c $BASE_DIR/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz) | cut -d " " -f 1 > $BASE_DIR/salmon_decoys.txt || { echo "Failed to create decoy list"; exit 1; }
sed -i.bak -e 's/>//g' $BASE_DIR/salmon_decoys.txt || { echo "Failed to modify decoy list"; exit 1; }
rm $BASE_DIR/salmon_decoys.txt.bak || { echo "Failed to remove backup decoy list"; exit 1; }

# Combine FASTA files
cat $BASE_DIR/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz $BASE_DIR/Rattus_norvegicus.mRatBN7.2.ncrna.fa.gz $BASE_DIR/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz > $BASE_DIR/Rattus_norvegicus.mRatBN7.2.gentrome.gz.fa || { echo "Failed to combine FASTA files"; exit 1; }

# Create Salmon index
salmon index -t $BASE_DIR/Rattus_norvegicus.mRatBN7.2.gentrome.gz.fa -d $BASE_DIR/salmon_decoys.txt -p 11 -i $BASE_DIR/mRatBN7_index || { echo "Failed to create Salmon index"; exit 1; }