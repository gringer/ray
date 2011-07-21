CODE=../code

g++ -I. -I$CODE $CODE/format/ColorSpaceCodec.cpp $CODE/core/OperatingSystem.cpp $CODE/memory/malloc_types.cpp $CODE/core/common_functions.cpp test_invert.cpp -g $CODE/cryptography/crypto.cpp -DMAXKMERLENGTH=32 $CODE/structures/Kmer.cpp -I ..
./a.out
