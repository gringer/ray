mpirun -np 30 ~/git-clones/ray/code/Ray -write-kmers -show-memory-usage \
-p /data/users/boiseb01/sra/solidsoftwaretools.com/ecoli2x50/ecoli_600x_F3.csfasta \
   /data/users/boiseb01/sra/solidsoftwaretools.com/ecoli2x50/ecoli_600x_R3.csfasta \
-o color-space-test
