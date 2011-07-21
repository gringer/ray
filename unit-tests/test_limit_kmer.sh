CODE=../code
g++ -g $CODE/format/ColorSpaceCodec.cpp $CODE/structures/Direction.cpp $CODE/core/OperatingSystem.cpp \
$CODE/structures/ReadAnnotation.cpp $CODE/memory/malloc_types.cpp $CODE/structures/Vertex.cpp test_kmer.cpp $CODE/core/common_functions.cpp $CODE/structures/Kmer.cpp $CODE/cryptography/crypto.cpp  -I. -I$CODE -D MAXKMERLENGTH=50 -DASSERT -I..
./a.out \
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
