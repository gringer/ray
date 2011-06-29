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
	int startPos = 0;
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


Read::Read(uint8_t*seq,int length){ // for raw sequences
	m_type=TYPE_SINGLE_END;
	m_sequence=seq;
	m_length=length;
	m_colorSpace = false;
	m_firstBase = 'N';
}

Read::Read(uint8_t*seq,int length,bool color){ // for raw sequences (colour space specified)
	m_type=TYPE_SINGLE_END;
	m_sequence=seq;
	m_length=length;
	m_colorSpace = color;
	m_firstBase = 'N';
}

Read::Read(uint8_t*seq,int length,char firstBase){ // for raw sequences with first base
	m_type=TYPE_SINGLE_END;
	m_sequence=seq;
	m_length=length;
	m_colorSpace = true;
	m_firstBase = firstBase;
}

Read::Read(const char*sequenceIn,MyAllocator*seqMyAllocator,bool trimFlag){
	string tSeq(sequenceIn);
	m_firstBase = CSC::bsChrToBS(tSeq.at(0));
	if(CSC::isColorSpace(tSeq)){
		m_colorSpace = true;
		// assume first character of sequence is first base, ensure this is [ACGTN]
		tSeq = tSeq.substr(1,tSeq.length()-1); // trim off first base
	} else {
		m_colorSpace = false;
	}
	m_forwardOffset=0;
	m_reverseOffset=0;
	m_type=TYPE_SINGLE_END;
	// clean sequence (convert to [ACGTN]+)
	tSeq = clean(tSeq);
	// trim sequence (remove initial / trailing 'N's)
	if(trimFlag){
		tSeq = trim(tSeq);
	}
	int length=tSeq.length();
	m_length=length;

	int requiredBytes=getRequiredBytes();

	uint8_t workingBuffer[requiredBytes];
	for(int i=0;i<requiredBytes;i++){
		workingBuffer[i]=0;
	}

	for(int position=0;position<length;position++){
		char nucleotide=tSeq.at(position);
		uint8_t code;
		if(m_colorSpace){
			code=CSC::csChrToInt(nucleotide);
		} else {
			code=CSC::bsChrToInt(nucleotide);
		}
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
		m_sequence=(uint8_t*)seqMyAllocator->allocate(requiredBytes*sizeof(uint8_t));
		memcpy(m_sequence,workingBuffer,requiredBytes);
	}
}

void Read::getSeq(char*workingBuffer,bool color,bool doubleEncoding) const{
	//TODO: get base-space sequence as colour-space
	if(!color && doubleEncoding){
		cout << "warning: useless double-encoding requested for base-space output... ";
	}
	string out("");
	out.reserve(m_length);
	for(int position=0;position<m_length;position++){
		int positionInWorkingBuffer=position/4;
		uint8_t word=m_sequence[positionInWorkingBuffer];
		int codePositionInWord=position%4;
		uint8_t code=(word<<(6-codePositionInWord*2));//eliminate bits before
		code=(code>>6);
		if(m_colorSpace){
			out += CSC::csIntToCS(code, doubleEncoding);
		} else {
			out += CSC::bsIntToBS(code);
		}
	}
	if(!color & m_colorSpace) {
		// output in base-space is desired, but read is in colour-space
		out.insert(0,1,m_firstBase);
		out = CSC::decodeCStoBS(out);
	} else {
		// output in colour-space is desired, but read is in base-space
		if(color & !m_colorSpace) {
			out = CSC::encodeBStoCS(out);
			out.erase(0,1); // remove first base to be consistent with current code
		}
		// TODO: append first base, once it it confirmed as 'safe' to do
		//out.insert(0,1,m_firstBase);
	}
	// Assumes workingBuffer has size 4000; this is consistent with the state as at 2011-06-28
	// 3999 allows \0 to be stored as final character
	strncpy(workingBuffer,out.c_str(),3999); // strncpy(dest, src, n)
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
	char buffer[4000];
	getSeq(buffer,color,false);
	return kmerAtPosition(buffer,pos,w,strand,color);
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

int Read::getRequiredBytes(){
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
	string   csConv1bs("TAGGGATATCTTTTTAATGCCGAAATAAGGGCTGATAACCGAGATTATTTC");
	string   csConv1cs("32002333220000303130320033020032123301032223033002");
	string csConv1csde("TGAAGTTTGGAAAATATCTATGAATTAGAATGCGTTACATGGGTATTAAG");
	Read colorSpaceRead2("T32002333220000303130320033..0032123301.3222303....", &ma, true);
	string   csConv2bs("TAGGGATATCTTTTTAATGCCGAAATAAAAATCAGCGGTTAGAGCCG");
	string   csConv2cs("3200233322000030313032003300003212330103222303");
	Read baseSpaceRead1("ACACCACGCAAAATATTTGCTCCAGCTCCTTTCATT", &ma, true);
	string   bsConv1cs("11101133100033300132201232202002130");
	Read baseSpaceRead2("ZZAC6ACGCXAAATAT55ACTCCAGCTCC..RCA..", &ma, true);
	// note: after cleaning and trimming, reads are stored internally,
	// replacing 'N' with 'A' or '0'.
	string   bsConv2cs("1101331000333300122012322010011");
	string   bsConv2bs("ACAACGCAAAATATAAACTCCAGCTCCAAACA");
	char buffer[4000];
	bool lastResult = true;
	bool allResults = true;
	cout << "1: checking colour-space encoding converted to double-encoded base-space... ";
	colorSpaceRead1.getSeq(buffer,false,true); // this may produce a warning
	lastResult = (csConv1bs.compare(buffer) == 0); allResults = allResults && lastResult;
	if(!lastResult) {
		cout << "Conversion failed\n" << buffer << endl
				<< csConv1bs << " expected\n";
	} else {
		cout << "success!\n";
	}
	cout << "2: checking colour-space encoding converted to colour-space... ";
	colorSpaceRead1.getSeq(buffer,true,false);
	lastResult = (csConv1cs.compare(buffer) == 0); allResults = allResults && lastResult;
	if(!lastResult) {
		cout << "Conversion failed\n" << buffer << endl
				<< csConv1cs << " expected\n";
	} else {
		cout << "success!\n";
	}
	cout << "3: checking colour-space encoding converted to double-encoded colour-space... ";
	colorSpaceRead1.getSeq(buffer,true,true);
	lastResult = (csConv1csde.compare(buffer) == 0); allResults = allResults && lastResult;
	if(!lastResult) {
		cout << "Conversion failed\n" << buffer << endl
				<< csConv1csde << " expected\n";
	} else {
		cout << "success!\n";
	}
	cout << "4: checking colour-space encoding with misreads converted to base-space... ";
	colorSpaceRead2.getSeq(buffer,false,false);
	lastResult = (csConv2bs.compare(buffer) == 0); allResults = allResults && lastResult;
	if(!lastResult) {
		cout << "Conversion failed\n" << buffer << endl
				<< csConv2bs << " expected\n";
	} else {
		cout << "success!\n";
	}
	cout << "5: checking colour-space encoding with misreads converted to colour-space... ";
	colorSpaceRead2.getSeq(buffer,true,false);
	lastResult = (csConv2cs.compare(buffer) == 0); allResults = allResults && lastResult;
	if(!lastResult) {
		cout << "Conversion failed\n" << buffer << endl
				<< csConv2cs << " expected\n";
	} else {
		cout << "success!\n";
	}
	cout << "6: checking base-space encoding converted to colour-space... ";
	baseSpaceRead1.getSeq(buffer,true,false);
	lastResult = (bsConv1cs.compare(buffer) == 0); allResults = allResults && lastResult;
	if(!lastResult) {
		cout << "Conversion failed\n" << buffer << endl
				<< bsConv1cs << " expected\n";
	} else {
		cout << "success!\n";
	}
	cout << "7: checking base-space encoding with misreads converted to colour-space... ";
	baseSpaceRead2.getSeq(buffer,true,false);
	lastResult = (bsConv2cs.compare(buffer) == 0); allResults = allResults && lastResult;
	if(!lastResult) {
		cout << "Conversion failed\n" << buffer << endl
				<< bsConv2cs << " expected\n";
	} else {
		cout << "success!\n";
	}
	cout << "7: checking base-space encoding with misreads converted to base-space... ";
	baseSpaceRead2.getSeq(buffer,false,false);
	lastResult = (bsConv2bs.compare(buffer) == 0); allResults = allResults && lastResult;
	if(!lastResult) {
		cout << "Conversion failed\n" << buffer << endl
				<< bsConv2bs << " expected\n";
	} else {
		cout << "success!\n";
	}
	return allResults;
}
