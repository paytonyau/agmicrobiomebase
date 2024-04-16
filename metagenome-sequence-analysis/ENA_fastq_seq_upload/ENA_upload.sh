#!/bin/bash
#SBATCH -o slurm-%x_%A.out   
#SBATCH --job-name=UpENA
#SBATCH --partition=long


#check for the correct number of args
if [ $# -ne 3 ]
then
    echo "Error in $0 - Invalid Argument Count"
    echo "Syntax: $0  <username> <password file> <file with paths to FASTQ files for transfer>"
    exit
fi

# Uploading files to ENA - European Nucleotide Archive

##
echo "starting file upload to ENA"
date

#ENA account username
username=$1

#file containing the ENA account password
passwd=`cat $2`

#array with the files for upload
#full paths, one per line
#read this in from an external config file supplied as a command line argument
files=`cat $3`

for file in $files
do
	echo "transferring file $file"
	date
	curl \
	--retry 12 \
	--retry-connrefused \
	-T $file \
	ftp://$username:$passwd@webin2.ebi.ac.uk
done

echo "file upload complete"
date


