CODE=../code

g++ -I. $CODE/format/ColorSpaceCodec.cpp $CODE/memory/malloc_types.cpp -I$CODE $CODE/core/common_functions.cpp test_invert.cpp -g $CODE/cryptography/crypto.cpp -DMAXKMERLENGTH=32 $CODE/structures/Kmer.cpp -I ..
./a.out
