CODE=../code

g++ $CODE/format/ColorSpaceCodec.cpp $CODE/structures/Direction.cpp $CODE/structures/ReadAnnotation.cpp $CODE/memory/malloc_types.cpp  $CODE/structures/Vertex.cpp test_kmer.cpp $CODE/core/OperatingSystem.cpp $CODE/core/common_functions.cpp $CODE/structures/Kmer.cpp $CODE/cryptography/crypto.cpp -I$CODE -I. -D MAXKMERLENGTH=32 -DASSERT -I..
./a.out TGAAATGGAAATGGTCTGGGAAG

g++ $CODE/format/ColorSpaceCodec.cpp $CODE/structures/Direction.cpp  $CODE/structures/ReadAnnotation.cpp $CODE/memory/malloc_types.cpp  $CODE/structures/Vertex.cpp test_kmer.cpp $CODE/core/OperatingSystem.cpp $CODE/core/common_functions.cpp $CODE/structures/Kmer.cpp $CODE/cryptography/crypto.cpp -I$CODE -I. -D MAXKMERLENGTH=64 -DASSERT -I..
./a.out \
TGAAATGGAAATGGTCTGGGAAAAACAACTAAAAGATATTATTGTAGTA

g++ $CODE/format/ColorSpaceCodec.cpp $CODE/structures/Direction.cpp  $CODE/structures/ReadAnnotation.cpp $CODE/memory/malloc_types.cpp  $CODE/structures/Vertex.cpp test_kmer.cpp $CODE/core/common_functions.cpp $CODE/core/OperatingSystem.cpp $CODE/structures/Kmer.cpp $CODE/cryptography/crypto.cpp -I$CODE -I. -D MAXKMERLENGTH=96 -DASSERT -I..
./a.out \
TGAAATGGAAATGGTCTGGGAAAAACAACTAAAAGATATTATTGTAGTAGCTGGTTTTGAAATTTATGACGCTGAAATAACTCCCCACTA

g++ $CODE/format/ColorSpaceCodec.cpp $CODE/structures/Direction.cpp  $CODE/structures/ReadAnnotation.cpp $CODE/memory/malloc_types.cpp  $CODE/structures/Vertex.cpp test_kmer.cpp $CODE/core/common_functions.cpp $CODE/structures/Kmer.cpp $CODE/core/OperatingSystem.cpp $CODE/cryptography/crypto.cpp -I$CODE -I. -D MAXKMERLENGTH=128 -DASSERT -I..
./a.out \
TGAAATGGAAATGGTCTGGGAAAAACAACTAAAAGATATTATTGTAGTAGCTGGTTTTGAAATTTATGACGCTGAAATAACTCCCCACTATATTTTCACCAAATTTATT

g++ $CODE/format/ColorSpaceCodec.cpp $CODE/structures/Direction.cpp  $CODE/structures/ReadAnnotation.cpp $CODE/memory/malloc_types.cpp  $CODE/structures/Vertex.cpp test_kmer.cpp $CODE/core/common_functions.cpp $CODE/structures/Kmer.cpp $CODE/cryptography/crypto.cpp $CODE/core/OperatingSystem.cpp -I$CODE -I. -D MAXKMERLENGTH=75 -DASSERT -I..
./a.out \
TGAAATGGAAATGGTCTGGGAAAAACAACTAAAAGATATTATTGTAGTAGCTGGTTTTGAAATTTATGACGCT

g++ $CODE/format/ColorSpaceCodec.cpp $CODE/structures/Direction.cpp  $CODE/structures/ReadAnnotation.cpp $CODE/memory/malloc_types.cpp  $CODE/structures/Vertex.cpp test_kmer.cpp $CODE/core/common_functions.cpp $CODE/structures/Kmer.cpp $CODE/cryptography/crypto.cpp $CODE/core/OperatingSystem.cpp -I$CODE -I. -D MAXKMERLENGTH=50 -DASSERT -I..
./a.out \
TGAAATGGAAATGGTCTGGGAAAAACAACTAAAAGATATTAT

CODE=../code; \
g++ $CODE/format/ColorSpaceCodec.cpp \
    $CODE/structures/Read.cpp \
    $CODE/core/common_functions.cpp \
    $CODE/memory/malloc_types.cpp \
    $CODE/memory/allocator.cpp \
    $CODE/memory/MyAllocator.cpp \
    $CODE/memory/ReusableMemoryStore.cpp \
    $CODE/structures/Kmer.cpp \
    $CODE/cryptography/crypto.cpp \
    $CODE/core/UnitTestHarness.cpp \
    unit_tests.cpp -I$CODE -I..
./a.out
