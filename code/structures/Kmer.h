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

#ifndef _Kmer
#define _Kmer

#include <core/constants.h>
#include <format/ColorSpaceCodec.h>
#include <stdint.h>
#include <vector>
#include <cstring>
#include <string>
#include <iostream>
#ifdef ASSERT
#include <assert.h>
#endif
#include <string>
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
#define KMER_STARTPIECE           (1)
#define KMER_DIRECTION_MASK       (0b1)
#define KMER_FORWARD_DIRECTION    (0b0)
#define KMER_REVERSE_DIRECTION    (0b1)
#define KMER_VALID                (0b10)
#define KMER_VALID_MASK           (0b10)
#define KMER_FLAGS_MASK           (0b11)
#define KMER_2BIT_MASK            (0b11)
#define KMER_CLEAR_FLAGS          (~((uint64_t)0b11))
#define KMER_LASTU64_LOCATION     (KMER_U64_ARRAY_SIZE - 1)

/**
 * Class for storing k-mers.
 * For now only an array of uint64_t is present, but later,
 * when the code is stable, this could be a mix of u64, u32 and u16 and u8
 * while maintaining the same interface, that are the two functions.
 *
 * k-mer format
 * bit 0: direction (0: forward, 1: reverse) -- always 0 for now
 * bit 1: validity (1: valid k-mer, 0: invalid)
 * bit 2 onwards: remaining sequence, 2 bits per base (in base-space)
 */

class Kmer{
	typedef ColorSpaceCodec CSC;
	/** the actual array of uint64_t */
	uint64_t m_u64[KMER_U64_ARRAY_SIZE];
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
	int compare(uint64_t* a, uint64_t* b)const;
	int compare(const Kmer& a)const;
	void printPieces();
	void pack(uint64_t*messageBuffer,int*messagePosition);
	void unpack(uint64_t*messageBuffer,int*messagePosition);
	void unpack(vector<uint64_t>*messageBuffer,int*messagePosition);
	string toBSString(int wordsize);
	/*
	 * transform a Kmer in a string
	 */
	string toString(int wordsize, bool showBases);
	Kmer rComp(int wordSize);
	uint8_t getFirstCode(int wordSize, bool asColorSpace);
	char getFirstSymbol(int wordSize, bool asColorSpace);
	uint8_t getLastCode(int wordSize, bool asColorSpace);
	char getLastSymbol(int wordSize, bool asColorSpace);
	int vertexRank(int arraySize,int wordSize);

	vector<Kmer> getEdges(uint8_t edges,int wordSize, bool outGoing);

	/**
	 * get the outgoing Kmer objects for a Kmer a having edges and
	 * a k-mer length k
	 */
	vector<Kmer> getOutgoingEdges(uint8_t edges,int wordSize);
	/**
	 * get the ingoing Kmer objects for a Kmer a having edges and
	 * a k-mer length k
	 */
	vector<Kmer> getIngoingEdges(uint8_t edges,int wordSize);
	Kmer& operator=(const Kmer&b);
	bool operator<(const Kmer&b)const;
	bool operator!=(const Kmer&b)const;
	bool operator==(const Kmer&b)const;

	INLINE
	void clear(){
		for(int i=0; i < KMER_U64_ARRAY_SIZE; i++){
			m_u64[i] = 0;
		}
	}

	INLINE
	void setPiece(int piece, int code, uint64_t* bitArray = NULL){
		if(bitArray == NULL){
			bitArray = m_u64;
		}
		int arrayPos = piece / 32;
		int bitLocation = (piece % 32) << 1;
		#ifdef ASSERT
		assert(code<4);
		assert(arrayPos < KMER_U64_ARRAY_SIZE);
		#endif
		// add in some range checking
		if(arrayPos < KMER_U64_ARRAY_SIZE){
			// note: cast to uint64_t is necessary to make shift work correctly
			bitArray[arrayPos] &= ~((uint64_t)KMER_2BIT_MASK << bitLocation); // clear previous set bits
			bitArray[arrayPos] |= ((uint64_t)code << bitLocation); // set new bits
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
			int code = ((m_u64[arrayPos] >> bitLocation) & KMER_2BIT_MASK);
			return(code); // retrieve bits from position
		} else {
			//cout << "error: attempt to access out of array bounds" << endl;
			return -1;
		}
	}

	INLINE
	int getPiece(int piece, uint64_t* bitArray){
		int arrayPos = piece / 32;
		int bitLocation = (piece % 32) << 1;
		#ifdef ASSERT
		assert(arrayPos < KMER_U64_ARRAY_SIZE);
		#endif
		if(arrayPos < KMER_U64_ARRAY_SIZE){
			int code = ((bitArray[arrayPos] >> bitLocation) & KMER_2BIT_MASK);
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

	/** hash 1 is used to distribute k-mers on MPI ranks */
	uint64_t hash_function_1();

	/** hash 2 is used for double hashing in the hash tables */
	uint64_t hash_function_2();

	bool isRC() const;

}ATTRIBUTE_PACKED;


#endif
