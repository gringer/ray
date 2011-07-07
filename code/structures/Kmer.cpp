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
#include <core/common_functions.h>
#include <format/ColorSpaceCodec.h>
#include <cryptography/crypto.h>

/*
 * initialise Kmer to a subsequence of a given string.
 */
Kmer::Kmer(string inSequence, int pos, int wordSize, char strand){
	// Based on previous code, if a sequence can fit in the number of
	// allocated uint64_t, then the transfer is allowed, even if the
	// sequence is greater than MAXKMERLENGTH
	// so, limit range checking to that done in set/get Piece

	// insert unknown first base, if necessary
	bool reverse = (strand == 'R')?true:false;
	bool validSequence = true;
	bool isCS = CSC::isColorSpace(inSequence);
	// note: all sequences are stored internally as colour-space
	// When storing colour-space, it is necessary to find the first base for the sub-sequence
	char lastBase = (inSequence.length()>0)?inSequence.at(0):'N';
	// an unknown first base, or no first base
	bool firstBaseKnown = (!isCS || (CSC::bsChrToBS(lastBase) != 'N'));
	/* if the input string (not the subsequence) starts in the middle of a
	 * colour-space sequence (e.g. '00120231'), then it is necessary to insert a 'N' */
	if(!firstBaseKnown && (CSC::csChrToInt(lastBase) > 3)){
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
	} else if(firstBaseKnown) {
		setPiece(1, CSC::bsChrToInt(lastBase));
	}
	int lastBaseCode = CSC::bsChrToInt(lastBase);
	int checkSum = 0;
	for(int i = 1; i < wordSize; i++){
		int nextCode = CSC::csChrToInt(inSequence.at(pos + i));
		if(!isCS){ // map that base to colour-space
			nextCode = CSC::mapBStoCS(lastBaseCode,nextCode);
		}
		if(nextCode > 3){ // uh oh... bad read
			cout << "Fatal: bad sequence: position= " << pos+i << endl;
			nextCode = 0;
			validSequence = false; // invalidate sequence
		}
		// get next base-space code (for reverse-complement, if necessary)
		lastBaseCode = CSC::mapCStoBS(lastBaseCode,nextCode);
		if(!reverse){
			setPiece(i+1, nextCode);
		} else {
			setPiece((wordSize - (i - 1)), nextCode);
		}
		checkSum = (checkSum + nextCode) % 4;
	}
	if(reverse){
		setPiece(1, lastBaseCode);
	}
	if(!validSequence){
		cout << "Fatal: invalid input sequence: " << inSequence << endl;
		firstBaseKnown = false; // make sure a checksum is used
		checkSum++; // sabotage the checksum
	}
	int flags = (KMER_COLORSPACE) +
			(firstBaseKnown?KMER_FIRSTBASE_KNOWN:0);
	setPiece(0,flags);
	if(!firstBaseKnown){
		setPiece(1,(4 - checkSum) % 4); // so a 2-bit checksum == 0
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
		assert(checkSum());
		#endif
}

Kmer::Kmer(const Kmer& b){
	clear();
	if(((b.m_u64[0] & (uint64_t)KMER_COLORSPACE) != 0)){
		// already in colour-space
		// so do a straight copy of data
		for(int i=0;i<getNumberOfU64();i++){
			m_u64[i] = b.m_u64[i];
		}
		#ifdef ASSERT
		assert(checkSum());
		#endif
	} else {
		// not in colour-space, so convert base-space to colour space
		int baseX = b.getPiece(1);
		setPiece(0,KMER_CS_FIRSTBASE_KNOWN);
		setPiece(1,baseX);
		for(int i=1;i<(MAXKMERLENGTH);i++){
			baseX = b.getPiece(i);
			int code = CSC::mapBStoCS(baseX,b.getPiece(i+1));
			code = (code>3)?0:code;
			setPiece(i+1,code);
		}
		// no checksum because first base is known
	}
}

Kmer::Kmer(){
	clear();
	// note: this Kmer has an invalid checksum, because it is not in
	// colour space, but has an unknown first base [00 00 .. .. ..]
}

Kmer::~Kmer(){
}

bool Kmer::checkSum(){
	// warn if first base is unknown and checksum is invalid
	// [this should alert to failed Kmer construction]
	if((m_u64[0] & KMER_FLAGS) == KMER_CS_FIRSTBASE_UNKNOWN){
		int checkSum = getPiece(1);
		for(int i=0;i<(MAXKMERLENGTH-1);i++){
			//+2: skip over flags and checksum / firstBase
			checkSum += getPiece(i+2);
		}
		return(checkSum == 0);
	} else if ((m_u64[0] & KMER_FLAGS) == KMER_BS_FIRSTBASE_UNKNOWN){
		// in base-space and first base unknown... not possible
		return false;
	} else {
		// no checkSum present, or not in colour-space, so assume the sequence is fine
		return true;
	}
}



int Kmer::getNumberOfU64(){
	return KMER_U64_ARRAY_SIZE;
}

void Kmer::printPieces(){
	for(int i=0;i<MAXKMERLENGTH;i++){
		printf("(%d)",getPiece(i));
	}
	printf("\n");
}

/*
 * edges is a 4-bit packed boolean array indicating which bases (or colours)
 * are ingoing edges for the Vertex that includes this Kmer. e.g. if A,G,T
 * (or 0,2,3) are all possible incoming edges (prefixes) for the vertex,
 * then the bit-packed representation is 0b1101
 *
 * This function generates a vector of all Kmers generated from inserting
 * bases/colours into the start of this Kmer (pushing a base off the end).
 *
 * If (for whatever reason) the first base is unknown in colour-space, this will
 * add '0' as the first colour, set the first base to unknown (rather than
 * the given base), and invalidate the Kmer checksum.
 */
vector<Kmer> Kmer::getIngoingEdges(uint8_t edges,int wordSize){
	vector<Kmer> inEdges;
	Kmer kmerTemplate;
	kmerTemplate.setPiece(0,getPiece(0)); // this copies the flags across
	int tempCheckSum = 0;
	bool isCS = isColorSpace();
	bool fbKnown = firstBaseKnown();
	int firstBase = getPiece(1);
	int startPiece = isCS?2:1;
	for(int i = startPiece; i < (wordSize); i++){
		// shift sequence 'up' one base/colour
		int nextPiece = getPiece(i);
		tempCheckSum = (tempCheckSum + nextPiece) % 4;
		kmerTemplate.setPiece(i+1, nextPiece);
	}
	for(int code = 0; code < 4; code++, edges >>= 1){
		if((edges & 1) == 1){
			// cout << "adding incoming edge " << CSC::bsIntToBS(code)
			// 	<< " to " << this->toString(wordSize, true)
			// 	<< " (" << this->toBSString(wordSize) << ") -> ";
			// copy from the template
			Kmer kmerIn(kmerTemplate);
			if(fbKnown){
				// add the new edge to the start of the Kmer
				kmerIn.setPiece(1,code);
				// recalculate the first colour-space code
			}
			if(isCS){
				int newCode = CSC::mapBStoCS(firstBase, code);
				if(newCode <= 3){
					kmerIn.setPiece(2,newCode);
				} else {
					kmerIn.setPiece(2,0);
					fbKnown = false;
					kmerIn.setPiece(2,tempCheckSum + 1); // sabotage checksum
				}
			}
			// cout << kmerIn.toString(wordSize, true) << endl;
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
 * This function is very similar to getIngoing edges:
 *   * sequence is shifted from (i+1) to (i), rather than (i) to (i+1)
 *   * outgoing edges are stored in upper 4 bits (so >> 4 prior to adding edges)
 *   * replacement piece is set is at wordSize, rather than sequence start
 *   * first base is recalculated by stepping one base forward
 *   * last base is necessary in order to calculate next colour (for colour-space)
 */
vector<Kmer> Kmer::getOutgoingEdges(uint8_t edges,int wordSize){
	vector<Kmer> outEdges;
	Kmer kmerTemplate;
	kmerTemplate.setPiece(0,getPiece(0)); // this copies the flags across
	int lastBaseCode = this->getLastCode(wordSize, false);
	int tempCheckSum = 0;
	bool isCS = isColorSpace();
	bool fbKnown = firstBaseKnown();
	int startPiece = isCS?2:1;
	for(int i = startPiece; i < (wordSize); i++){
		// shift sequence 'down' one base/colour
		int nextPiece = getPiece(i+1);
		tempCheckSum = (tempCheckSum + nextPiece) % 4;
		kmerTemplate.setPiece(i, nextPiece);
	}
	// outgoing edges are stored in the high bits, so shift down 4
	edges >>= 4;
	for(int code = 0; code < 4; code++, edges >>= 1){
		int newCheckSum = tempCheckSum;
		if((edges & 1) == 1){
			// cout << "adding outgoing edge " << CSC::bsIntToBS(code)
			// 	<< " to " << this->toString(wordSize, true)
			// 	<< " (" << this->toBSString(wordSize) << ") -> ";
			// copy from the template
			Kmer kmerOut(kmerTemplate);
			// add the new edge to the end of the Kmer
			int newCode = code;
			if(isCS){
				newCode = CSC::mapBStoCS(lastBaseCode,code);
				if(newCode > 3){
					newCode = 0;
					fbKnown = false;
					newCheckSum++; // sabotage the checksum
				}
			}
			kmerOut.setPiece(wordSize,newCode);
			// place first base / checksum (for colour-space Kmers)
			if(isCS){
				if(!fbKnown){
					// first base unknown, so generate checksum
					kmerOut.setPiece(1,(4 - ((newCheckSum + newCode) % 4)) % 4);
				} else {
					// recalculate first base by stepping forward 1
					// A0X -> [0]|AX, A1X -> [1]|CX, etc.
					kmerOut.setPiece(1, CSC::mapCStoBS(getPiece(1),getPiece(2)));
				}
			}
			// cout << kmerOut.toString(wordSize, true) << endl;
			// finally, send to vector
			outEdges.push_back(kmerOut);
		}
	}
	return outEdges;
}


int Kmer::getFirstCode(bool asColorSpace) {
	#ifdef ASSERT
	assert(checkSum());
	#endif
	bool isCS = isColorSpace();
	int firstCode = 4;
	if(asColorSpace){
			if(isCS){
				// piece 0 flags, piece 1 firstBase...
				// colour-space sequence begins from piece 2
				firstCode = getPiece(2);
			} else {
				firstCode = CSC::mapBStoCS(getPiece(1),getPiece(2));
			}
	} else {
		bool firstBaseKnown = ((getPiece(0) & KMER_FIRSTBASE_KNOWN) != 0);
		if(firstBaseKnown){
			firstCode = getPiece(1);
		} else {
			// first base unknown, so [0-3] are all equally bad
			// and piece 1 contains the checksum
			firstCode = 4; // well... you did ask...
		}
	}
	return(firstCode);
}


char Kmer::getFirstSymbol(bool asColorSpace){
	if(asColorSpace){
		bool firstBaseKnown = ((getPiece(0) & KMER_FIRSTBASE_KNOWN) != 0);
		if(!firstBaseKnown){
			return 'N';
		}
		return CSC::csIntToCS(getFirstCode(asColorSpace),false);
	} else {
		return CSC::bsIntToBS(getFirstCode(asColorSpace));
	}
}

int Kmer::getLastCode(int wordSize, bool asColorSpace){
	#ifdef ASSERT
	assert(checkSum());
	#endif
	bool isCS = isColorSpace();
	if(asColorSpace == isCS){
		if(isCS){
			// add 1 to skip over first base
			return(getPiece(wordSize+1));
		} else {
			// sequence starts at piece 1
			return(getPiece(wordSize));
		}
	}
	bool firstBaseKnown = ((getPiece(0) & KMER_FIRSTBASE_KNOWN) != 0);
	if(asColorSpace && !firstBaseKnown){
		return(4);
	}
	// these will be somewhat expensive to calculate
	int lastBaseCode = getPiece(1);
	for(int i = 1; i < (wordSize); i++){
		int code = getPiece(i+1);
		if(isCS){
			lastBaseCode = CSC::mapCStoBS(lastBaseCode,code);
		} else {
			lastBaseCode = CSC::mapBStoCS(lastBaseCode,code);
		}
	}
	return lastBaseCode;
}


char Kmer::getLastSymbol(int wordSize, bool asColorSpace){
	if(asColorSpace){
		if(!this->firstBaseKnown()){
			return 'N';
		}
		return CSC::csIntToCS(getLastCode(wordSize,asColorSpace),false);
	} else {
		return CSC::bsIntToBS(getLastCode(wordSize,asColorSpace));
	}
}

string Kmer::toString(int wordSize, bool showFirstBase){
	// Based on previous code, if the asked-for length fits in the number of
	// allocated uint64_t, then the request is allowed, even if the
	// sequence is greater than MAXKMERLENGTH (so things beyond that length
	// should be undefined)
	// so, limit range checking to that done in set/get Piece
	#ifdef ASSERT
	assert(checkSum());
	#endif
	string out("");
	out.reserve(wordSize);
	bool isCS = isColorSpace();
	bool fbKnown = firstBaseKnown();
	for(int i = 0; i < wordSize; i++){
		int code = getPiece(i+1);
		if(i == 0){
			if(showFirstBase || !isCS){
				// in base space, want to show everything
				if(!fbKnown){
					out += "N";
				} else {
					out += CSC::bsIntToBS(code);
				}
			}
		} else {
			if(!isCS){
				out += CSC::bsIntToBS(code);
			} else {
				out += CSC::csIntToCS(code, false);
			}
		}
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
	assert(checkSum());
	assert(wordSize <= KMER_MAX_PIECES);
	#endif
	string out("");
	out.reserve(wordSize);
	bool isCS = isColorSpace();
	bool fbKnown = firstBaseKnown();
	if(!fbKnown){
		out.append(wordSize,'N'); // first base unknown, so save a few processing steps
	} else {
		int lastBase = getPiece(1);
		out += CSC::bsIntToBS(lastBase);
		for(int i = 1; i < (wordSize); i++){
			if(!isCS){
				out += CSC::bsIntToBS(getPiece(i+1));
			} else {
				int code = getPiece(i+1);
				lastBase = CSC::mapCStoBS(lastBase,code);
				out += CSC::bsIntToBS(lastBase);
			}
		}
	}
	return out;
}

/*
 * Return the reverse-complement of this Kmer. The reverse complement will
 * be returned as a colour-space Kmer.
 */

Kmer Kmer::rComp(int wordSize){
	// Based on previous code, if a sequence can fit in the number of
	// allocated uint64_t, then the transfer is allowed, even if the
	// sequence is greater than MAXKMERLENGTH
	// so, limit range checking to that done in set/get Piece
	#ifdef ASSERT
	assert(checkSum());
	#endif
	int flags = getPiece(0);
	bool colorSpace = ((flags & KMER_COLORSPACE) != 0);
	bool firstBaseKnown = ((flags & KMER_FIRSTBASE_KNOWN) != 0);
	if(!colorSpace){
		Kmer tK(*this); // make a new copy as Colour space
		return(tK.rComp(wordSize));
	}
	Kmer tRC;
	int lastBase = firstBaseKnown?CSC::complement(getPiece(1)):4;
	int checkSum = 0;
	// it is necessary to track through and convert to base
	// space in order to get the reverse-complement first base
	for(int piece = 0; piece < (wordSize-1); piece++){
		int code = getPiece(piece+2);
		checkSum = (checkSum + code) % 4;
		lastBase = CSC::mapCStoBS(lastBase,code);
		tRC.setPiece(wordSize - piece, code);
	}
	if(lastBase > 3){
		tRC.setPiece(0, KMER_CS_FIRSTBASE_UNKNOWN);
		tRC.setPiece(1, (4 - checkSum) % 4);
	} else {
		tRC.setPiece(0,KMER_CS_FIRSTBASE_KNOWN);
		tRC.setPiece(1, lastBase);
	}
	return tRC;
}

void Kmer::pack(uint64_t*messageBuffer,int*messagePosition){
	for(int i=0;i<getNumberOfU64();i++){
		messageBuffer[*messagePosition]=getU64(i);
		(*messagePosition)++;
	}
}

void Kmer::unpack(uint64_t*messageBuffer,int*messagePosition){
	for(int i=0;i<getNumberOfU64();i++){
		setU64(i,messageBuffer[*messagePosition]);
		(*messagePosition)++;
	}
}

void Kmer::unpack(vector<uint64_t>*messageBuffer,int*messagePosition){
	for(int i=0;i<getNumberOfU64();i++){
		setU64(i,(*messageBuffer)[*messagePosition]);
		(*messagePosition)++;
	}
}

uint64_t Kmer::getHash_1(){
	uint64_t key=m_u64[0];
	for(int i=0;i<getNumberOfU64();i++){
		uint64_t newKey = m_u64[i];
		if(i == 0){
			// clear unknown bit and first base, as these may change
			newKey &= KMER_CLEAR_FIRSTBASE;
		}
		key ^= uniform_hashing_function_1_64_64(newKey);
	}
	return key;
}

uint64_t Kmer::getHash_2(){
	uint64_t key=m_u64[0];
	for(int i=0;i<getNumberOfU64();i++){
		uint64_t newKey = m_u64[i];
		if(i == 0){
			// clear unknown bit and first base, as these may change
			newKey &= KMER_CLEAR_FIRSTBASE;
		}
		key ^= uniform_hashing_function_2_64_64(newKey);
	}
	return key;
}

void Kmer::operator=(const Kmer&b){
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		m_u64[i]=b.m_u64[i];
	}
}

bool Kmer::isColorSpace() const{
	return((m_u64[0] & (uint64_t)KMER_COLORSPACE) != 0);
}

bool Kmer::firstBaseKnown() const{
	return((m_u64[0] & (uint64_t)KMER_FIRSTBASE_KNOWN) != 0);
}


int Kmer::compare(const Kmer& b) const{
	// normalise to make sure both Kmers are in the same space
	if(this->isColorSpace() && !b.isColorSpace()){
		Kmer csB(b);
		return(this->compare(csB));
	}
	if(!this->isColorSpace() && b.isColorSpace()){
		Kmer csA(*this);
		return(csA.compare(b));
	}
	// compare sequence at all known locations
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		uint64_t checkA = m_u64[i];
		uint64_t checkB = b.m_u64[i];
		if(i == 0){
			if(((checkA & KMER_FIRSTBASE_KNOWN) == 0) ||
					((checkB & KMER_FIRSTBASE_KNOWN) == 0)){
				// at least one of the first bases is unknown, so ignore first base for comparison
				checkA &= KMER_CLEAR_FIRSTBASE; // clear unknown bit and first base
				checkB &= KMER_CLEAR_FIRSTBASE;
			}
		}
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
	return (this->compare(b) < 0);
}

bool Kmer::isLower(Kmer*a){
	return (this->compare(*a) < 0);
}

bool Kmer::isEqual(Kmer*a){
	return (this->compare(*a) == 0);
}

