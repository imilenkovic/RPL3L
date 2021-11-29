## HMMsearch

# Download the hmm profile from the pfam domain database

# Download the proteomes of species of interest (Uniprot -> proteomes)

# hmmsearch for the profile within the proteome

#!/bin/bash
for i in /Users/ivan/Bioinfo/hmm/Proteomes2/*fasta;
do
hmmsearch -o ${i}.hmmsearchoutput Ribosomal_L3.hmm $i
done

# Move all the .hmmsearchoutput files to a new folder
# Extract the highest hits from the output files

#!/bin/bash
for i in /Users/ivan/Bioinfo/hmm/29072019_RPL3/*output;
do
head -23 $i | tail -6 > $i.ID
done

# Concatenate all files into a single one
cat *.ID > trial1
# Get rid of mitochondrial proteins (the word can be cut in the output file so it's enough to use the first few letters)
sed '/Mitoch/d' trial1 > trial2
sed '/mitoch/d' trial2 > trial3
# Get rid of fragmented protein inputs
sed '/Fragment/d' trial3 > trial4
# Get rid of redundant lines
sed '/>/d' trial4 > trial5
sed '/Domain/d' trial5 > trial6
sed '/inclusion/d' trial6 > trial7
sed '/^$/d' trial7 > trial8
# Extract the IDs and gene names
cut -c 59-83 trial8 > trial9
# Replace | with space to facilitate copying into excel
sed 's/|/ /g' trial9 > trial10
# Copy the trial10 table into excel and paste special by space as a delimiter
# Make three columns, "Uniprot" (the ID), "GeneName" (from the trial10 table) and "Species" (manually add)
# Make new directory to perform hmmalign
mkdir FastaUniprot
cd FastaUniprot
# Write the excel table into a file (with vim or nano) with the headers
vim Uniprot_GeneName_species
# Write another file with only the protein IDs (no header)
vim uniprot.list
# The script to replace Uniprot ID with Gene name and Species

use strict;
use warnings;
my %mem=();
open A, $ARGV[0];
while (<A>) {
	chomp;
	my @a = split /\s+/,$_;
	$mem{$a[0]} = $a[1]."_".$a[2];
}

open B, $ARGV[1];

while (<B>) {
	if (/^>/) {
		my @a = split /\|/,$_;
		$a[1] = $mem{$a[1]} if exists $mem{$a[1]};
		my $l = join "|",@a;
		print ">$a[1]\n" ;
	}
	else {
		print $_;
	}
}

# Script to extract fasta from uniprot
#%%%%%%%%Extracting Fasta from Uniprot lists of families%%%%%%%
#extract_fasta.sh
#Extracting fasta sequences for the proteins of interest in a list of orthologs proteins of a gene group
#Dependencies: perl
#Notes: replace.pl (added) is a perl script to replace Uniprot ID with Gene and species name by using a file called Uniprot_Species_file 
#Uniprot_Species_file contains a table of 3 columns with Uniprot ID, GeneName and Species

#Create new directory where to work
mkdir dir.$1
cd dir.$1
cp ../$1 .

# List of IDs to uniprot fasta
cat $1 | while read line; do
     #wget http://www.uniprot.org/uniprot/$line.fasta
     curl -L -O http://www.uniprot.org/uniprot/$line.fasta > $line.fasta
done

# Merge all fasta sequences
cat *fasta > $1.fasta

#change the name of fasta to GeneName_Species (replace.pl and Uniprot_GeneName_Species files are in the folder)
perl ../replace.pl ../Uniprot_GeneName_Species $1.fasta > $1.named.fasta

#Run by typing:
extract_fasta.sh uniprotlist.file

# Do hmmalign with the uniprot.list.named.fasta file (copy the hmm profile into the folder):
hmmalign Ribosomal_L3.hmm uniprot.list.named.fasta > 29072019_align

# Convert the stockholm format to fasta format 
http://sequenceconversion.bugaco.com/converter/biology/sequences/stockholm_to_fasta.php


# Build the tree with iqtree 

iqtree -s 29072019_align.fasta -st AA -m TEST -bb 5000 -alrt 5000 -nt 8 -minsup 0.5

# Submit the job
qsub -cwd -pe smp 8 -l virtual_free=4G rpl3l.sh
