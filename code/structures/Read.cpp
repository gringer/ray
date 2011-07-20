/*
 	Ray
    Copyright (C)  2010, 2011  Sébastien Boisvert

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

#include<assert.h>
#include<core/common_functions.h>
#include<structures/Read.h>
#include<cstdlib>
#include<iostream>
#include<cstring>
#include<format/ColorSpaceCodec.h>
#include<core/UnitTestHarness.h>
using namespace  std;

/*
 * Convert a sequence to [ACGTN]+. Colour-space characters are double-encoded as tr/0123/ACGT/.
 * All other characters are converted to 'N'.
 */
string Read::clean(string sequence){
	string out(sequence);
	//cout<<"   In='"<<out<<"'"<<endl;
	for(string::iterator it = out.begin(); it < out.end(); it++){
		// convert to [ACGTN]
		if(m_colorSpace){
			*it = CSC::csChrToDE(*it);
		} else {
			*it = CSC::bsChrToBS(*it);
		}
	}
	//cout<<"Clean='"<<out<<"'"<<endl;
	return out;
}

/*
 * Trim off 'N' from the start and end of the sequence.
 */
string Read::trim(string sequence){
	string out(sequence);
	//cout<<"In='"<<out<<"'"<<endl;
	// discard N at the beginning and end of the read.
	// erase up to the first symbol that is a A,T,C or G
	uint8_t startPos = 0;
	while((startPos < out.length()) && (out.at(startPos) == 'N')){
		startPos++;
	}
	out.erase(0,startPos);
	//cout<<"Trimmed first "<<out<<endl;
	// erase from the string end back to the last symbol that is a A,T,C, or G
	int endPos = out.length();
	while((endPos > 0) && (out.at(endPos-1) == 'N')){
		endPos--;
	}
	out.erase(endPos, out.length() - endPos);
	//cout<<"Trimmed last "<<out<<endl;
	return out;
}

Read::Read(){ // needed for repeated reads (assembler/seedExtender.cpp:1034)
	m_type=TYPE_SINGLE_END;
	m_length=0;
}


Read::Read(uint8_t*seq,int length,bool color, bool firstBaseKnown){ // for raw sequences
	m_type=TYPE_SINGLE_END;
	m_sequence=seq;
	m_length=length;
	m_colorSpace = color;
	m_firstBaseKnown = firstBaseKnown;
}

Read::Read(const Read& b,MyAllocator*seqMyAllocator){
	// note: trimming the read doesn't make sense, as the internal representation for b is 2-bit
	m_length = b.m_length;
	m_type = b.m_type;

	m_colorSpace = b.m_colorSpace;
	m_firstBaseKnown = b.m_firstBaseKnown;

	m_forwardOffset = b.m_forwardOffset;
	m_reverseOffset = b.m_reverseOffset;

	m_pairedRead = b.m_pairedRead; //TODO: is it safe to do this?

	// copy over the sequence data, allocating new bytes as necessary
	int requiredBytes = getRequiredBytes();
	if(b.m_sequence==NULL){
		m_sequence=NULL;
	}else{
		if(seqMyAllocator == NULL){ // used to simplify unit tests
			m_sequence = new uint8_t[requiredBytes];
		} else {
			m_sequence = (uint8_t*)seqMyAllocator->allocate(requiredBytes*sizeof(uint8_t));
		}
		memcpy(m_sequence,b.m_sequence,requiredBytes);
	}
}


Read::Read(string sequenceIn,MyAllocator*seqMyAllocator,bool trimFlag){
	char firstBase = CSC::bsChrToBS(sequenceIn.at(0));
	bool inColorSpace = CSC::isColorSpace(sequenceIn);
	if(inColorSpace){
		if(CSC::csChrToInt(firstBase) != CSC::bsChrToInt(firstBase)){
			// first character is a colour-space code, so insert unknown first base
			sequenceIn.insert(0,1,'N');
			firstBase = 'N';
		}
	}
	// all sequences are stored internally as colorSpace
	m_firstBaseKnown = (firstBase != 'N');
	m_colorSpace = true;
	m_forwardOffset=0;
	m_reverseOffset=0;
	m_type=TYPE_SINGLE_END;
	// clean sequence (convert to [ACGTN]+)
	sequenceIn = clean(sequenceIn);
	// trim sequence (remove initial / trailing 'N's)
	if(trimFlag){
		sequenceIn = trim(sequenceIn);
		if(!inColorSpace && sequenceIn.length()>0){
			firstBase = sequenceIn.at(0);
			m_firstBaseKnown = true;
		} else if(!m_firstBaseKnown){
			sequenceIn.insert(0,1,'A'); // insert dummy first base because it was trimmed off
		}
	}
	m_length=sequenceIn.length();
	if(!inColorSpace){
		sequenceIn = CSC::encodeBStoCS(sequenceIn);
	}

	int requiredBytes=getRequiredBytes();

	uint8_t workingBuffer[requiredBytes];
	for(int i=0;i<requiredBytes;i++){
		workingBuffer[i]=0;
	}
	for(int position=0;position<m_length;position++){
		char nucleotide=sequenceIn.at(position);
		uint8_t code=CSC::csChrToInt(nucleotide);
		if(code > 3){
			code=0; // TODO: note that this introduces error bias. It would be better to
		}	        // split reads into contiguous sequences of unambiguous bases
		#ifdef __READ_VERBOSITY
		if(position%4==0){
			cout<<"|";
		}
		cout<<" "<<(int)code;
		#endif
		int positionInWorkingBuffer=position/4;
		int codePositionInWord=position%4;
		uint8_t wordToUpdate=workingBuffer[positionInWorkingBuffer];
		// shift the code and or with the word to update
		code=(code<<(codePositionInWord*2));
		wordToUpdate=wordToUpdate|code;
		workingBuffer[positionInWorkingBuffer]=wordToUpdate;
	}
	#ifdef __READ_VERBOSITY
	cout<<endl;
	for(int i=0;i<requiredBytes;i++){
		cout<<" "<<(int)workingBuffer[i];
	}

	cout<<endl;
	#endif
	if(requiredBytes==0){
		m_sequence=NULL;
	}else{
		if(seqMyAllocator == NULL){ // used to simplify unit tests
			m_sequence = new uint8_t[requiredBytes];
		} else {
			m_sequence=(uint8_t*)seqMyAllocator->allocate(requiredBytes*sizeof(uint8_t));
		}
		memcpy(m_sequence,workingBuffer,requiredBytes);
	}
}

string Read::getSeq(bool color,bool doubleEncoding) const{
	if(!color && doubleEncoding){
		cout << "warning: useless double-encoding requested for base-space output... ";
	}
	string out("");
	out.reserve(m_length);
	for(int position=0;position<m_length;position++){
		int positionInWorkingBuffer=position/4;
		uint8_t word=m_sequence[positionInWorkingBuffer];
		int bitPositionInWord=(position%4) * 2;
		uint8_t code=(word >> bitPositionInWord) & 0b11; // mask out word
		if(position == 0){
			out += m_firstBaseKnown?CSC::bsIntToBS(code):'N';
		} else {
			out += CSC::csIntToCS(code, doubleEncoding);
		}
	}
	if(!color) {
		// output in base-space is desired, but read is in colour-space
		out = CSC::decodeCStoBS(out);
	}
	return out;
}

int Read::length()const{
	return m_length;
}

/*                      
 *           -----------------------------------
 *           -----------------------------------
 *                     p p-1 p-2               0
 */
Kmer Read::getVertex(int pos,int w,char strand,bool color) const {
	//TODO: possibly do a copy of bytes, without the intermediate string
	return Kmer(getSeq(color,false),pos,w,strand);
}

bool Read::hasPairedRead()const{
	return m_type!=TYPE_SINGLE_END;
}

PairedRead*Read::getPairedRead(){
	if(m_type==TYPE_SINGLE_END){
		return NULL;
	}
	return &m_pairedRead;
}

uint8_t*Read::getRawSequence(){
	return m_sequence;
}

int Read::getRequiredBytes() const{
	int requiredBits=2*m_length;
	int modulo=requiredBits%8;
	if(modulo!=0){
		int bitsToAdd=8-modulo;
		requiredBits+=bitsToAdd;
	}

	#ifdef ASSERT
	assert(requiredBits%8==0);
	#endif

	int requiredBytes=requiredBits/8;
	return requiredBytes;
}
void Read::setLeftType(){
	m_type=TYPE_LEFT_END;
}

void Read::setRightType(){
	m_type=TYPE_RIGHT_END;
}

int Read::getType(){
	return m_type;
}

void Read::setType(uint8_t a){
	m_type=a;
}

void Read::setForwardOffset(int a){
	m_forwardOffset=a;
}

void Read::setReverseOffset(int a){
	m_reverseOffset=a;
}

int Read::getForwardOffset(){
	return m_forwardOffset;
}

int Read::getReverseOffset(){
	return m_reverseOffset;
}

bool Read::check(){
	MyAllocator ma;
	ma.constructor(4194304,14,false); // from Loader::constructor
	Read colorSpaceRead1("T32002333220000303130320033020032123301032223033002", &ma, true);
	Read colorSpaceRead2("T32002333220000303130320033..0032123301.3222303....", &ma, true);
	Read baseSpaceRead1("ACACCACGCAAAATATTTGCTCCAGCTCCTTTCATT", &ma, true);
	Read baseSpaceRead2("ZZAC6ACGCXAAATAT55ACTCCAGCTCC..RCA..", &ma, true);
	// note: after cleaning and trimming, reads are stored internally,
	// replacing 'N' with 'A' or '0'.
	UnitTestHarness uth("Read");
	string sequence;
	uth.preProcessTest("colour-space encoding converted to double-encoded base-space");
	sequence = colorSpaceRead1.getSeq(false,true); // this may produce a warning
	uth.compareOutput(sequence, "TAGGGATATCTTTTTAATGCCGAAATAAGGGCTGATAACCGAGATTATTTC");
	uth.preProcessTest("colour-space encoding converted to colour-space");
	sequence = colorSpaceRead1.getSeq(true,false);
	uth.compareOutput(sequence, "T32002333220000303130320033020032123301032223033002");
	uth.preProcessTest("colour-space encoding converted to double-encoded colour-space");
	sequence = colorSpaceRead1.getSeq(true,true);
	uth.compareOutput(sequence, "TTGAAGTTTGGAAAATATCTATGAATTAGAATGCGTTACATGGGTATTAAG");
	uth.preProcessTest("colour-space encoding with misreads converted to base-space");
	sequence = colorSpaceRead2.getSeq(false,false);
	uth.compareOutput(sequence, "TAGGGATATCTTTTTAATGCCGAAATAAAAATCAGCGGTTAGAGCCG");
	uth.preProcessTest("colour-space encoding with misreads converted to colour-space");
	sequence = colorSpaceRead2.getSeq(true,false);
	uth.compareOutput(sequence, "T3200233322000030313032003300003212330103222303");
	uth.preProcessTest("base-space encoding converted to colour-space");
	sequence = baseSpaceRead1.getSeq(true,false);
	uth.compareOutput(sequence, "A11101133100033300132201232202002130");
	uth.preProcessTest("base-space encoding with misreads converted to colour-space");
	sequence = baseSpaceRead2.getSeq(true,false);
	uth.compareOutput(sequence, "A1001330000333000122012322000001");
	uth.preProcessTest("base-space encoding with misreads converted to base-space");
	sequence = baseSpaceRead2.getSeq(false,false);
	uth.compareOutput(sequence, "ACCCATAAAAATATTTTGAGGTCGAGGGGGGT");
	return uth.getSuccess();
}
