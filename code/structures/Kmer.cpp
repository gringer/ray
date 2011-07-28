/*
 	Ray
    Copyright (C) 2011  Sébastien Boisvert

	http://DeNovoAssembler.SourceForge.Net/

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You have received a copy of the GNU General Public License
    along with this program (COPYING).  
	see <http://www.gnu.org/licenses/>

*/

#include <structures/Kmer.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <string>
#include <core/common_functions.h>
#include <format/ColorSpaceCodec.h>
#include <cryptography/crypto.h>

using namespace std;

/*
 * initialise Kmer to a subsequence of a given string.
 */
Kmer::Kmer(string inSequence, int pos, int wordSize, char strand){
	// Based on previous code, if a sequence can fit in the number of
	// allocated uint64_t, then the transfer is allowed, even if the
	// sequence is greater than MAXKMERLENGTH
	// so, limit range checking to that done in set/get Piece

	// TODO: store only the lower of forward/RC, and change strand/direction flag appropriately

	// insert unknown first base, if necessary
	bool reverse = (strand == 'R')?true:false;
	bool validSequence = true;
	/*
	 * * all sequences are stored internally as base-space
	 * * sequence.at(1) makes sure range checking happens
	 * * it is assumed that any colour-space sequence will
	 *   have a colour at its second position
	 */
	bool isCS = ((CSC::csChrToInt(inSequence.at(1)) != CSC::bsChrToInt(inSequence.at(1))));

	// When loading from colour-space, it is necessary to find the first base for the sub-sequence
	char lastBase = (inSequence.length()>0)?inSequence.at(0):'N';
	// an unknown first base, or no first base
	bool firstBaseKnown = (!isCS || (CSC::bsChrToBS(lastBase) != 'N'));
	/* if the input string (not the subsequence) starts in the middle of a
	 * colour-space sequence (e.g. '00120231'), then it is necessary to insert a 'N' */
	if(!firstBaseKnown && (CSC::csChrToInt(lastBase) <= 3)){
		inSequence.insert(0,1,'N');
		lastBase = 'N';
	}

	// check and adjust ranges
	if(!reverse && (strand != 'F')){
		cout << "Fatal: strand direction is incorrect: strand= " << strand << endl;
		validSequence = false;
		pos = 0;
		wordSize = 0; // invalid strand, don't write any more sequence
	}
	if(pos >= (int)inSequence.length()){
		cout << "Fatal: offset is too large: position= " << pos << " Length= "
				<< inSequence.length() << " WordSize=" << wordSize << endl;
		validSequence = false;
		pos = 0;
		wordSize = 0;
	}
	if((pos + wordSize) > inSequence.length()){
		// don't read in more than the available sequence
		wordSize = inSequence.length() - pos;
	}
	if(wordSize < 0){
		cout << "Fatal: invalid word size: position= " << pos << " Length= "
				<< inSequence.length() << " WordSize=" << wordSize << endl;
		validSequence = false;
		wordSize = 0;
	}
	if(wordSize > KMER_MAX_PIECES){
		cout << "Fatal: word size cannot fit in Kmer array: position= " << pos << " Length= "
				<< inSequence.length() << " WordSize=" << wordSize << endl;
		validSequence = false;
		wordSize = 0;
	}
	#ifdef ASSERT
	assert(validSequence);
	#endif
	// finished with the range checking, and variable setup. Now onto conversion
	clear();

	//note: inSequence will be at least 1 character long, because 'N' is inserted for an empty string
	string priorDecode = isCS?CSC::decodeCStoBS(inSequence.substr(0,pos+1), false):inSequence.substr(pos,1);
	lastBase = priorDecode.at(priorDecode.length()-1);
	if(CSC::bsChrToBS(lastBase) == 'N'){
		firstBaseKnown = false;
	}
	if(reverse){
		lastBase = complementNucleotide(toupper(lastBase));
	}
	int lastBaseCode = CSC::bsChrToInt(lastBase);
	for(int i = 1; i < wordSize; i++){
		int nextCode = CSC::csChrToInt(inSequence.at(pos + i));
		if(!isCS){ // map that base to colour-space
			nextCode = CSC::mapBStoCS(lastBaseCode,nextCode);
		}
		if(lastBaseCode > 3){ // uh oh... bad read
			cout << "Fatal: bad sequence: position= " << pos+i << endl;
			lastBaseCode = 0;
			validSequence = false; // invalidate sequence
		}
//		cout << "Last base code: " << CSC::bsIntToBS(lastBaseCode) << endl;
		if(!reverse){
			setPiece(i, lastBaseCode);
		} else {
			setPiece((wordSize - (i - 1)), lastBaseCode);
		}
		// get next base-space code
		lastBaseCode = CSC::mapCStoBS(lastBaseCode,nextCode);
	}
	// add in final base
	if(lastBaseCode > 3){
		lastBaseCode = 0;
		validSequence = false; // invalidate sequence
	} else {
		setPiece(wordSize, lastBaseCode);
	}
	if(!validSequence){
		cout << "Fatal: invalid input sequence: " << inSequence << endl;
		#ifdef ASSERT
		assert(validSequence);
		#endif
		clear();
	} else {
		int flags = KMER_FORWARD_DIRECTION | (validSequence?KMER_VALID:0);
		setPiece(0,flags);
	}
}

/*
 * Used for initializing from raw bit sequences (e.g. messages)
 */
Kmer::Kmer(uint64_t* rawBits){
		for(int i=0;i<getNumberOfU64();i++){
			m_u64[i]=rawBits[i];
		}
		#ifdef ASSERT
		assert(isValid());
		#endif
}

Kmer::Kmer(const Kmer& b){
	for(int i=0;i<getNumberOfU64();i++){
		m_u64[i] = b.m_u64[i];
	}
	//TODO: no validity check because EdgePurgeWorker calls this for a brand new k-mer on array access
}

Kmer::Kmer(){
	clear();
	// note: this k-mer is invalid, because the valid flag is not set [00 .. .. ..]
}   //                                                                  *

Kmer::~Kmer(){
}

/* check if the k-mer is valid
 * [this should alert to failed k-mer construction]
 */
bool Kmer::isValid(){
	return((m_u64[0] & KMER_VALID_MASK) != 0);
}

int Kmer::getNumberOfU64(){
	return KMER_U64_ARRAY_SIZE;
}

void Kmer::printPieces() const {
	printf("[flags=%d]",getPiece(0));
	for(int i=KMER_STARTPIECE;i<MAXKMERLENGTH;i++){
		printf("(%d)",getPiece(i));
	}
	printf("\n");
}

/*
 * edges is a 4-bit packed boolean array indicating which bases (or colours)
 * are ingoing edges for the Vertex that includes this k-mer. e.g. if A,G,T
 * (or 0,2,3) are all possible incoming edges (prefixes) for the vertex,
 * then the bit-packed representation is 0b1101
 *
 * This function generates a vector of all k-mers generated from inserting
 * bases/colours into the start of this k-mer (pushing a base off the end).
 */
vector<Kmer> Kmer::getIngoingEdges(uint8_t edges,int wordSize){
	vector<Kmer> inEdges;
	Kmer kmerTemplate;
	kmerTemplate.setPiece(0,getPiece(0)); // this copies the flags across
	for(int i = KMER_STARTPIECE; i < (wordSize); i++){
		// shift sequence 'up' one base/colour
		int nextPiece = getPiece(i);
		kmerTemplate.setPiece(i+1, nextPiece);
	}
	for(int code = 0; code < 4; code++, edges >>= 1){
		if((edges & 1) == 1){
			// cout << "adding incoming edge " << CSC::bsIntToBS(code)
			// 	<< " to " << this->toString(wordSize, true)
			// 	<< " (" << this->toBSString(wordSize) << ") -> ";
			// copy from the template
			Kmer kmerIn(kmerTemplate);
			// add the new edge to the start of the Kmer
			kmerIn.setPiece(KMER_STARTPIECE,code);
			// finally, send to vector
			inEdges.push_back(kmerIn);
		}
	}
	return inEdges;
}

/*
 * edges is a 4-bit packed boolean array indicating which bases (or colours)
 * are ingoing edges for the Vertex that includes this Kmer. e.g. if A,G,T
 * (or 0,2,3) are all possible incoming edges (prefixes) for the vertex,
 * then the bit-packed representation is 0b1101
 *
 * This function generates a vector of all Kmers generated from inserting
 * bases/colours at the end of this Kmer (pushing a base off the start).
 *
 * This function is very similar to getIngoing edges, except for these points:
 *   * sequence is shifted from (i+1) to (i), rather than (i) to (i+1)
 *   * outgoing edges are stored in upper 4 bits (so >> 4 prior to adding edges)
 *   * replacement piece is set is at wordSize, rather than sequence start
 */
vector<Kmer> Kmer::getOutgoingEdges(uint8_t edges,int wordSize){
	vector<Kmer> outEdges;
	Kmer kmerTemplate;
	kmerTemplate.setPiece(0,getPiece(0)); // this copies the flags across
	for(int i = KMER_STARTPIECE; i < (wordSize); i++){
		// shift sequence 'down' one base/colour
		int nextPiece = getPiece(i+1);
		kmerTemplate.setPiece(i, nextPiece);
	}
	// outgoing edges are stored in the high bits, so shift down 4
	edges >>= 4;
	for(int code = 0; code < 4; code++, edges >>= 1){
		if((edges & 1) == 1){
			// cout << "adding outgoing edge " << CSC::bsIntToBS(code)
			// 	<< " to " << this->toString(wordSize, true)
			// 	<< " (" << this->toBSString(wordSize) << ") -> ";
			// copy from the template
			Kmer kmerOut(kmerTemplate);
			// add the new edge to the end of the Kmer
			kmerOut.setPiece(wordSize,code);
			// finally, send to vector
			outEdges.push_back(kmerOut);
		}
	}
	return outEdges;
}


uint8_t Kmer::getFirstCode(bool asColorSpace) {
	#ifdef ASSERT
	assert(isValid());
	#endif
	if(asColorSpace){
		return CSC::mapBStoCS(getPiece(KMER_STARTPIECE),getPiece(KMER_STARTPIECE+1));
	} else {
		return getPiece(KMER_STARTPIECE);
	}
}


char Kmer::getFirstSymbol(bool asColorSpace){
	if(asColorSpace){
		return CSC::csIntToCS(getFirstCode(asColorSpace));
	} else {
		return CSC::bsIntToBS(getFirstCode(asColorSpace));
	}
}

uint8_t Kmer::getLastCode(int wordSize, bool asColorSpace){
	#ifdef ASSERT
	assert(isValid());
	#endif
	if(asColorSpace){
		return CSC::mapBStoCS(getPiece(wordSize-1),getPiece(wordSize));
	} else {
		return getPiece(wordSize);
	}
}

char Kmer::getLastSymbol(int wordSize, bool asColorSpace){
	if(asColorSpace){
		return CSC::csIntToCS(getLastCode(wordSize,asColorSpace));
	} else {
		return CSC::bsIntToBS(getLastCode(wordSize,asColorSpace));
	}
}

string Kmer::toString(int wordSize, bool showBases){
	// Based on previous code, if the asked-for length fits in the number of
	// allocated uint64_t, then the request is allowed, even if the
	// sequence is greater than MAXKMERLENGTH (so things beyond that length
	// should be undefined)
	// so, limit range checking to that done in set/get Piece
	#ifdef ASSERT
	assert(isValid());
	#endif
	string out("");
	if(showBases){
		out += CSC::bsIntToBS(getPiece(KMER_STARTPIECE));
	}
	// output colour-space sequence
	int oldBase = getPiece(KMER_STARTPIECE);
	for(int i = KMER_STARTPIECE+1; i <= wordSize; i++){
		int nextBase = getPiece(i);
		out += CSC::csIntToCS(CSC::mapBStoCS(oldBase, nextBase));
		oldBase = nextBase;
	}
	if(showBases){
		out += CSC::bsIntToBS(getPiece(wordSize));
	}
	return out;
}

string Kmer::toBSString(int wordSize){
	if(wordSize > KMER_MAX_PIECES){
		cout << "Warning: wordSize (" << wordSize << ") " << "is greater than maximum allowed ("
				<< KMER_MAX_PIECES << "). Kmer output will be trimmed to maximum allowable size" << endl;
		wordSize = KMER_MAX_PIECES;
	}
	#ifdef ASSERT
	assert(isValid());
	assert(wordSize <= KMER_MAX_PIECES);
	#endif
	string out("");
	for(int i = KMER_STARTPIECE; i <= wordSize; i++){
		out += CSC::bsIntToBS(getPiece(i));
	}
	return out;
}

/*
 * Return the reverse-complement of this Kmer. The reverse complement will
 * be returned as a colour-space Kmer.
 * TODO: reverse complement should be constant time, or at most O(getNumberOfU64()), rather than O(wordSize)
 */

Kmer Kmer::rComp(int wordSize){
	// Based on previous code, if a sequence can fit in the number of
	// allocated uint64_t, then the transfer is allowed, even if the
	// sequence is greater than MAXKMERLENGTH
	// so, limit range checking to that done in set/get Piece
	#ifdef ASSERT
	assert(isValid());
	#endif
	int flags = getPiece(0);
	Kmer tRC;
	tRC.setPiece(0, flags);
	// copy pieces across in reverse order
	int copied = 0;
	for(int piece = 0; piece < wordSize; piece++){
		int code = getPiece(KMER_STARTPIECE + piece);
		tRC.setPiece(wordSize - piece, CSC::complement(code));
		copied++;
	}
	return tRC;
}

void Kmer::pack(uint64_t*messageBuffer,int*messagePosition){
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		messageBuffer[*messagePosition]=getRawBits(i);
		(*messagePosition)++;
	}
}

void Kmer::unpack(uint64_t*messageBuffer,int*messagePosition){
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		m_u64[i] = messageBuffer[*messagePosition];
		(*messagePosition)++;
	}
}

void Kmer::unpack(vector<uint64_t>*messageBuffer,int*messagePosition){
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		m_u64[i] = messageBuffer->at(*messagePosition);
		(*messagePosition)++;
	}
}

uint64_t Kmer::hash_function_1(){
	uint64_t key=m_u64[0];
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		uint64_t newKey = m_u64[i];
		if(i == 0){
			// clear flags, as these are not relevant to the sequence
			newKey &= KMER_CLEAR_FLAGS;
		}
		key ^= uniform_hashing_function_1_64_64(newKey);
	}
	return key;
}

uint64_t Kmer::hash_function_2(){
	uint64_t key=m_u64[0];
	for(int i=0;i<getNumberOfU64();i++){
		uint64_t newKey = m_u64[i];
		if(i == 0){
			// clear unknown bit and first base, as these may change
			newKey &= KMER_CLEAR_FLAGS;
		}
		key ^= uniform_hashing_function_2_64_64(newKey);
	}
	return key;
}

Kmer& Kmer::operator=(const Kmer&b){
	if(this != &b){
		for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
			m_u64[i]=b.m_u64[i];
		}
		#ifdef ASSERT
		if(!isValid()){
			cout << "Address: " << &b << endl;
			b.printPieces();
			this->printPieces();
			assert(isValid());
		}
		#endif
	}
	return *this;
}

int Kmer::compare(const Kmer& b) const{
	// compare sequence at all known locations
	// this moves in reverse so that the most significant base for comparison is the last base
	for(int i=KMER_U64_ARRAY_SIZE-1;i>=0;i--){
		uint64_t checkA = m_u64[i];
		uint64_t checkB = b.m_u64[i];
		if(checkA<checkB){
			return -1;
		}else if(checkA>checkB){
			return 1;
		}
	}
	return 0;
}

bool Kmer::operator==(const Kmer&b)const{
	return (this->compare(b) == 0);
}

bool Kmer::operator!=(const Kmer&b)const{
	return (this->compare(b) != 0);
}

bool Kmer::operator<(const Kmer&b)const{
	return (compare(b) < 0);
}

bool Kmer::isLower(Kmer*a){
	return (compare(*a) < 0);
}

bool Kmer::isEqual(Kmer*a){
	return (compare(*a) == 0);
}

int Kmer::vertexRank(int arraySize,int wordSize){
	Kmer b=rComp(wordSize);
	//TODO: work out why the copy is necessary here
	if(isLower(&b)){
		b=*this;
	}
	return b.hash_function_1()%(arraySize);
}
