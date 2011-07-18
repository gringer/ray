#!/bin/sh
echo -n "Checking full Ray run with phiX genome... ";
../readSimulator/VirtualNextGenSequencer tests/phix/phix.fasta 0 200 10 5000 50 tests_phix/1.fasta tests/phix/2.fasta > /dev/null && echo -n "1000 Reads simulated... "
echo -n "Running Ray... ";
mpirun -np 2  ../code/Ray --debug-seeds -p tests/phix/1.fasta tests/phix/2.fasta -o ray_output/test_phiX
fasta_formatter -i tests/phix/phix.fasta | grep $(fasta_formatter -i ray_output/test_phiX.Scaffolds.fasta | grep -v '^>') > /dev/null && echo "success (match in forward direction)!"
fasta_formatter -i tests/phix/phix.fasta | fastx_reverse_complement | grep $(fasta_formatter -i ray_output/test_phiX.Scaffolds.fasta | grep -v '^>') > /dev/null && echo "success (match in reverse direction)!"
