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
#include <stdio.h>
#include <assert.h>

Kmer::Kmer(string sequence){
	#ifdef ASSERT
	assert(sequence.length() == MAXKMERLENGTH);
	#endif
	clear();
	int checkSum = 0;
	bool colorSpace = CSC::isColorSpace(sequence);
	bool firstBaseKnown = (!colorSpace || (sequence.at(i) != 'N'));
	for(int i = 0; i < sequence.length(); i++){
		int code = charToCode(sequence.at(i));
		if(i > 0){ // only checksum the non-firstBase sequence
			checkSum = (checkSum + code) % 4;
		}
		setPiece(i+1, code);
	}
	int flags = (colorSpace?1:0) + (firstBaseKnown?2:0);
	setPiece(0,flags);
	if(!firstBaseknown){
		setPiece(1,(4 - checkSum) % 4); // so a 2-bit checksum == 0
	}
}

Kmer::Kmer(){
	clear();
}

Kmer::~Kmer(){
}

int Kmer::getNumberOfU64(){
	return KMER_U64_ARRAY_SIZE;
}

bool Kmer::isLower(Kmer*a){
	return (this < (*a));
}

bool Kmer::isEqual(Kmer*a){
	return (this == (*a));
}

void Kmer::print(){
	for(int i=0;i<getNumberOfU64();i++){
		printf("(%b)",m_u64[i]);
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

bool Kmer::operator<(const Kmer&b)const{
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		checkA = m_u64[i];
		checkB = b.m_u64[i];
		if(i == 0){
			if((checkA & 0b10 == 0) || (checkB & 0b10 == 0)){
				// at least one of the first bases is unknown, so ignore first base for comparison
				checkA &= ~(0b1110); // clear unknown bit and first base
				checkB &= ~(0b1110);
			}
		}
		if(checkA<checkB){
			return true;
		}else if(checkA>checkB){
			return false;
		}
	}
	return false;
}

void Kmer::operator=(const Kmer&b){
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		m_u64[i]=b.m_u64[i];
	}
}

bool Kmer::operator==(const Kmer&b) const{
	//TODO: make this work for base-space vs colour-space
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		checkA = m_u64[i];
		checkB = b.m_u64[i];
		if(i == 0){
			if((checkA & 0b10 == 0) || (checkB & 0b10 == 0)){
				// at least one of the first bases is unknown, so ignore first base for comparison
				checkA &= ~(0b1110); // clear unknown bit and first base
				checkB &= ~(0b1110);
			}
		}
		if(m_u64[i]!=b.m_u64[i]){
			return false;
		}
	}
	return true;
}

bool Kmer::operator!=(const Kmer&b) const{
	//TODO: make this work for base-space vs colour-space
	for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
		checkA = m_u64[i];
		checkB = b.m_u64[i];
		if(i == 0){
			if((checkA & 0b10 == 0) || (checkB & 0b10 == 0)){
				// at least one of the first bases is unknown, so ignore first base for comparison
				checkA &= ~(0b1110); // clear unknown bit and first base
				checkB &= ~(0b1110);
			}
		}
		if(m_u64[i]!=b.m_u64[i]){
			return true;
		}
	}
	return false;
}

