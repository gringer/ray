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
#ifdef ASSERT
#include <assert.h>
#endif
using namespace std;

/*
 * Determine the number of uint64_t necessary to host 
 * k-mers of maximum length MAXKMERLENGTH
 */
#define KMER_REQUIRED_BITS (2*MAXKMERLENGTH+2)
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

/**
 * Class for storing k-mers.
 * For now only an array of uint64_t is present, but later,
 * when the code is stable, this could be a mix of u64, u32 and u16 and u8
 * while maintening the same interface, that are the two functions.
 *
 */
class Kmer{
	typedef ColorSpaceCodec CSC;
	uint64_t m_u64[KMER_U64_ARRAY_SIZE];
	//TODO could possibly also store as a bitset
public:
	Kmer();
	~Kmer();
	bool isEqual(Kmer*a);
	bool isLower(Kmer*a);
	void print();
	void pack(uint64_t*messageBuffer,int*messagePosition);
	void unpack(uint64_t*messageBuffer,int*messagePosition);
	void unpack(vector<uint64_t>*messageBuffer,int*messagePosition);
	void operator=(const Kmer&b);
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
		m_u64[arrayPos] &= ~(0b11 << bitLocation); // clear previous set bits
		m_u64[arrayPos] |= (code << bitLocation); // set new bits
	}

	INLINE
	void setU64(int i,uint64_t b){
		#ifdef ASSERT
		assert(i<KMER_U64_ARRAY_SIZE);
		#endif
		m_u64[i]=b;
	}

	INLINE
	uint64_t getU64(int i){
		#ifdef ASSERT
		assert(i<KMER_U64_ARRAY_SIZE);
		#endif
		return m_u64[i];
	}

	int getNumberOfU64();

}ATTRIBUTE_PACKED;


#endif
