/*
 	Ray
    Copyright (C)  2010  Sébastien Boisvert

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

#ifndef _Read
#define _Read

#include <string>
#include <stdint.h>
#include <vector>
#include <format/ColorSpaceCodec.h>
#include <memory/MyAllocator.h>
#include <structures/PairedRead.h>
using namespace std;

#define TYPE_SINGLE_END 0
#define TYPE_LEFT_END 1
#define TYPE_RIGHT_END 2

/**
 * a read is represented as a uint8_t*,
 * 2 bits per nucleotide
 * and a (possible) link to paired information.
 *
 * Note: reads are all stored internally as colour-space
 */
class Read{
	uint8_t *m_sequence;
	PairedRead m_pairedRead;// the read on the left

	/* maximum value: 65535 */
	uint16_t m_length;
	uint8_t m_type;
	
	// for colour-space representation
	bool m_colorSpace; // should always be true
	bool m_firstBaseKnown; // first base is stored in first bit position

	// for the scaffolder:
	uint8_t m_forwardOffset;
	uint8_t m_reverseOffset;

	typedef ColorSpaceCodec CSC; // makes calling static functions a bit less painful

	string clean(string sequence);
	string trim(string sequence);
public:
	Read(); // needed for repeated reads (assembler/seedExtender.cpp:1034)
	// raw sequence
	Read(uint8_t*seq,int length,bool color = false, bool firstBaseKnown = true);
	Read(const Read& b,MyAllocator*seqMyAllocator);
	Read(string sequenceIn,MyAllocator*seqMyAllocator,bool trim);
	string getSeq(bool color)const;
	int length()const;
	Kmer getVertex(int pos,int w,char strand,bool color)const;
	bool hasPairedRead()const;
	PairedRead*getPairedRead();
	uint8_t*getRawSequence();
	int getRequiredBytes()const;
	void setRawSequence(uint8_t*seq,int length);
	void setRightType();
	void setLeftType();
	int getType();
	void setType(uint8_t type);
	void setForwardOffset(int a);
	void setReverseOffset(int a);
	int getForwardOffset();
	int getReverseOffset();
	static bool check();
} ATTRIBUTE_PACKED;

#endif
