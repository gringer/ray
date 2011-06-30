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
#include <core/common_functions.h>
#include <format/ColorSpaceCodec.h>

Kmer::Kmer(string sequence){
	#ifdef ASSERT
	// note: a sequence AGCT will be indistinguishable from AGCTAAAAAAAAAA in base-space
	//       and AGCTTTTTTTTTTTT in colour-space
	assert(sequence.length() <= MAXKMERLENGTH);
	#endif
	clear();
	int checkSum = 0;
	bool colorSpace = CSC::isColorSpace(sequence);
	bool firstBaseKnown = (!colorSpace || (sequence.at(0) != 'N'));
	for(int i = 0; i < (int)sequence.length(); i++){
		int code = CSC::csChrToInt(sequence.at(i));
		if(i > 0){ // only checksum the non-firstBase sequence
			checkSum = (checkSum + code) % 4;
		}
		setPiece(i+1, code>3?0:code);
	}
	int flags = (colorSpace?1:0) + (firstBaseKnown?2:0);
	setPiece(0,flags);
	if(!firstBaseKnown){
		setPiece(1,(4 - checkSum) % 4); // so a 2-bit checksum == 0
	}
}

Kmer::Kmer(const Kmer& b, bool convertToColourSpace){
	clear();
	if((b.m_u64[0] & 1) == 1){
		// already in colour-space, so do a straight copy of the data
		for(int i=0;i<getNumberOfU64();i++){
			m_u64[i] = b.m_u64[i];
		}
		checkSum();
	} else {
		// not in colour-space, so convert base-space to colour space
		int baseX = b.getPiece(1);
		setPiece(0,KMER_BS_FIRSTBASE_KNOWN);
		setPiece(1,baseX);
		for(int i=2;i<(MAXKMERLENGTH);i++){
			int code = CSC::mapBStoCS(baseX,b.getPiece(i));
			setPiece(i,code>3?0:code);
		}
	}
}

Kmer::Kmer(){
	clear();
}

Kmer::~Kmer(){
}

bool Kmer::checkSum(){
	// warn if first base is unknown and checksum is invalid
	// [this should alert to incorrect encodings... eventually]
	if((m_u64[0] & KMER_FLAGS) == KMER_CS_FIRSTBASE_UNKNOWN){
		int checkSum = (m_u64[0] & KMER_FIRSTBASE) >> 2;
		for(int i=0;i<(MAXKMERLENGTH-1);i++){
			//+2: skip over flags and checksum / firstBase
			checkSum += getPiece(i+2);
		}
		#ifdef ASSERT
		assert(checkSum == 0);
		#endif
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

void Kmer::print(){
	for(int i=0;i<MAXKMERLENGTH;i++){
		printf("(%d)",getPiece(i));
	}
	printf("\n");
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

void Kmer::operator=(const Kmer&b){
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		m_u64[i]=b.m_u64[i];
	}
}

bool Kmer::isColorSpace() const{
	return((m_u64[0] & 1) == 1);
}

int Kmer::compare(const Kmer& b) const{
	// normalise to make sure both Kmers are in the same space
	if(this->isColorSpace() && !b.isColorSpace()){
		Kmer csB(b, true);
		return(this->compare(csB));
	}
	if(!this->isColorSpace() && b.isColorSpace()){
		Kmer csA(*this, true);
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

