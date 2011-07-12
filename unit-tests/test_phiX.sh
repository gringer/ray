#!/bin/sh
echo -n "Checking full Ray run with phiX genome... ";
genomeID=$(wget -q 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&tool=ray&term=NC_001422' -O - | grep '<Id>' | perl -pe 's/\s*<.*?>//g')
wget -q "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&tool=ray&rettype=fasta&id="${genomeID} -O phix_genome.fasta && echo -n "phiX sequence downloaded... "
../readSimulator/VirtualNextGenSequencer phix_genome.fasta 0 200 10 1000 50 phix_1.fasta phix_2.fasta > /dev/null && echo -n "1000 Reads simulated... "
echo -n "Running Ray... ";
mpirun -np 2  ../code/Ray --debug-seeds -p phix_1.fasta phix_2.fasta
fasta_formatter -i phix_genome.fasta | grep $(fasta_formatter -i RayOutput.Scaffolds.fasta | grep -v '^>') > /dev/null && echo "success (match in forward direction)!"
fasta_formatter -i phix_genome.fasta | fastx_reverse_complement | grep $(fasta_formatter -i RayOutput.Scaffolds.fasta | grep -v '^>') > /dev/null && echo "success (match in reverse direction)!"
