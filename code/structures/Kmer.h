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

/* TODO: merge the k-mer-related things here */

#ifndef _Kmer
#define _Kmer

#include <core/constants.h>
#include <format/ColorSpaceCodec.h>
#include <stdint.h>
#include <vector>
#include <string>
#include <iostream>
#ifdef ASSERT
#include <assert.h>
#endif
using namespace std;

/*
 * Determine the number of uint64_t necessary to host 
 * k-mers of maximum length MAXKMERLENGTH
 */
#define KMER_REQUIRED_BITS (2*MAXKMERLENGTH+4)
#define KMER_REQUIRED_BYTES (KMER_REQUIRED_BITS/8)
#define KMER_REQUIRED_BYTES_MODULO (KMER_REQUIRED_BITS%8)
#if KMER_REQUIRED_BYTES_MODULO
	#define KMER_BYTES (KMER_REQUIRED_BYTES+1)
#else
	#define KMER_BYTES KMER_REQUIRED_BYTES
#endif

#define KMER_UINT64_T (KMER_BYTES/8)
#define KMER_UINT64_T_MODULO (KMER_BYTES%8)
#if KMER_UINT64_T_MODULO
	#define KMER_U64_ARRAY_SIZE (KMER_UINT64_T+1)
#else
	#define KMER_U64_ARRAY_SIZE (KMER_UINT64_T)
#endif

#define KMER_MAX_BITS (KMER_U64_ARRAY_SIZE*64)
#define KMER_MAX_PIECES (KMER_MAX_BITS/2 - 2)

//#define KMER_CS_FIRSTBASE_UNKNOWN (0b01)
//#define KMER_CS_FIRSTBASE_KNOWN   (0b11)
//#define KMER_BS_FIRSTBASE_UNKNOWN (0b00)
//#define KMER_BS_FIRSTBASE_KNOWN   (0b10)
//#define KMER_COLORSPACE           (0b01)
// piece 0 flags, piece 1 firstBase...
// colour-space sequence begins from piece 2
#define KMER_FIRSTBASE_LOCATION   (1)
#define KMER_LASTBASE_LOCATION    (KMER_MAX_PIECES + 1)
#define KMER_STARTPIECE           (2)
#define KMER_DIRECTION            (0b0)
#define KMER_FORWARD_DIRECTION    (0b0)
#define KMER_REVERSE_DIRECTION    (0b1)
#define KMER_FIRSTBASE_KNOWN      (0b10)
#define KMER_FLAGS                (0b11)
#define KMER_2BITMASK             (0b11)
#define KMER_FIRSTBASE            (0b1100)
#define KMER_FIRSTBASE_UNKNOWN_F  (0b11)
#define KMER_CLEAR_FIRSTBASE      (~(KMER_FIRSTBASE))
#define KMER_CLEAR_LASTBASE       (~((uint64_t)0b11 << 62))
#define KMER_LASTU64_LOCATION     (KMER_U64_ARRAY_SIZE - 1)

/**
 * Class for storing k-mers.
 * For now only an array of uint64_t is present, but later,
 * when the code is stable, this could be a mix of u64, u32 and u16 and u8
 * while maintaining the same interface, that are the two functions.
 *
 * Kmer format
 * bit 0: direction (0: forward, 1: reverse) -- always 0 for now
 * bit 1: first base known (1: known first base, 0: unknown first base)
 * bit 2-3: first base (A=0b00,C=0b01,G=0b10,T=0b11)
 *        Note: if first-base is unknown, this should be 3 (0b11)
 * bit 4 onwards: remaining sequence, 2 bits per base/colour
 * last 2 bits: last base (in colour-space)
 *        Note: if first-base is unknown, this should be 3 (0b11)
 *
 */

class Kmer{
	typedef ColorSpaceCodec CSC;
	uint64_t m_u64[KMER_U64_ARRAY_SIZE];
	//TODO could possibly also store this as a bitset
public:
	//note: all constructors convert to colour-space
	Kmer();
	Kmer(string inSequence, int pos = 0, int wordSize = KMER_MAX_PIECES, char strand = 'F');
	Kmer(uint64_t* rawBits);
	Kmer(const Kmer& b);
	~Kmer();
	bool isValid();
	bool isEqual(Kmer*a);
	bool isLower(Kmer*a);
	int compare(const Kmer& a)const;
	void printPieces();
	void pack(uint64_t*messageBuffer,int*messagePosition);
	void unpack(uint64_t*messageBuffer,int*messagePosition);
	void unpack(vector<uint64_t>*messageBuffer,int*messagePosition);
	uint64_t getHash_1();
	uint64_t getHash_2();
	string toBSString(int wordsize);
	string toString(int wordsize, bool showBases);
	Kmer rComp(int wordSize);
	vector<Kmer> getIngoingEdges(uint8_t edges,int wordSize);
	vector<Kmer> getOutgoingEdges(uint8_t edges,int wordSize);
	int getFirstCode(bool asColorSpace);
	char getFirstSymbol(bool asColorSpace);
	int getLastCode(int wordSize, bool asColorSpace);
	char getLastSymbol(int wordSize, bool asColorSpace);
	void operator=(const Kmer&b);
	bool isColorSpace() const;
	bool firstBaseKnown() const;
	bool operator<(const Kmer&b)const;
	bool operator!=(const Kmer&b)const;
	bool operator==(const Kmer&b)const;

	INLINE
	void clear(){
		for(int i=0;i<getNumberOfU64();i++){
			m_u64[i]=0;
		}
	}

	INLINE
	void setPiece(int piece, int code){
		int arrayPos = piece / 32;
		int bitLocation = (piece % 32) << 1;
		#ifdef ASSERT
		assert(code<4);
		assert(arrayPos < KMER_U64_ARRAY_SIZE);
		#endif
		// add in some range checking
		if(arrayPos < KMER_U64_ARRAY_SIZE){
			// note: cast to uint64_t is necessary to make shift work correctly
			m_u64[arrayPos] &= ~((uint64_t)KMER_2BITMASK << bitLocation); // clear previous set bits
			m_u64[arrayPos] |= ((uint64_t)code << bitLocation); // set new bits
		}
	}

	INLINE
	int getPiece(int piece) const{
		int arrayPos = piece / 32;
		int bitLocation = (piece % 32) << 1;
		#ifdef ASSERT
		assert(arrayPos < KMER_U64_ARRAY_SIZE);
		#endif
		if(arrayPos < KMER_U64_ARRAY_SIZE){
			int code = ((m_u64[arrayPos] >> bitLocation) & KMER_2BITMASK);
			return(code); // retrieve bits from position
		} else {
			//cout << "error: attempt to access out of array bounds" << endl;
			return -1;
		}
	}

	int getNumberOfU64();
	INLINE
	uint64_t getRawBits(int arrayLocation) const{
		#ifdef ASSERT
		assert(arrayLocation<KMER_U64_ARRAY_SIZE);
		#endif
		return m_u64[arrayLocation];
	}


private:
	INLINE
	void setU64(int i,uint64_t b){
		#ifdef ASSERT
		assert(i<KMER_U64_ARRAY_SIZE);
		#endif
		m_u64[i]=b;
	}

	INLINE
	uint64_t getU64(int i) const{
		#ifdef ASSERT
		assert(i<KMER_U64_ARRAY_SIZE);
		#endif
		return m_u64[i];
	}

}ATTRIBUTE_PACKED;


#endif
