/*
 	Ray
    Copyright (C) 2010, 2011  Sébastien Boisvert

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

#include <cryptography/crypto.h>
#include <assert.h>
#include <stdio.h>
#include <core/constants.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <core/common_functions.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>

#ifdef OS_POSIX
#include <unistd.h>
#endif
#ifdef OS_WIN
#include <windows.h>
#endif

using namespace std;

char complementNucleotide(char c){
	switch(c){
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
		default:
			return c;
	}
}

string reverseComplement(string*a){
	ostringstream b;
	for(int i=a->length()-1;i>=0;i--){
		b<<complementNucleotide((*a)[i]);
	}
	return b.str();
}

bool isValidDNA(string x){
	for(string::iterator i=x.begin();i<x.end();i++){
		 // colour-space is valid DNA as well
		if(CSC::csChrToInt(*i) > 3){
			return false;
		}
	}
	return true;
}

string addLineBreaks(string dna,int columns){
	ostringstream output;
	int j=0;
	while(j<(int)dna.length()){
		output<<dna.substr(j,columns)<<endl;
		j+=columns;
	}
	return output.str();
}

string convertToString(vector<Kmer>*b,int m_wordSize,bool color){
	ostringstream a;
	#ifdef USE_DISTANT_SEGMENTS_GRAPH
	//
	//TTAATT
	// TTAATT
	//  TTAATT
	//  the first vertex can not fill in the first delta letter alone, it needs help.
	for(int p=0;p<m_wordSize;p++){
		a<<(b->at(p)).getFirstSymbol(color);
	}
	#else
	 //TODO: output in base-space for now, because there may be scaffolder problems otherwise
	a<<((*b)[0]).toBSString(m_wordSize);
	#endif
	for(int j=1;j<(int)(*b).size();j++){
		a<<(*b)[j].getLastSymbol(m_wordSize,false);
	}
	string contig=a.str();
	return contig;
}

int vertexRank(Kmer*a,int _size,int w,bool color){
	Kmer b=a->rComp(w);
	if(b.isLower(a)){
		return(b.hash_function_1() % _size);
	} else {
		return(a->hash_function_1() % _size);
	}
}

int roundNumber(int s,int alignment){
	return ((s/alignment)+1)*alignment;
}

/** real-time */
uint64_t getMilliSeconds(){
	uint64_t milliSeconds=0;
	#ifdef HAVE_CLOCK_GETTIME
	timespec temp;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&temp);
	uint64_t seconds=temp.tv_sec;
	uint64_t nanoseconds=temp.tv_nsec;
	milliSeconds=seconds*1000+nanoseconds/1000/1000;
	#endif
	return milliSeconds;
}

void showMemoryUsage(int rank){
	#ifdef __linux__
	ifstream f("/proc/self/status");
	while(!f.eof()){
		string key;
		f>>key;
		if(key=="VmData:"){
			uint64_t count;
			f>>count;
			cout<<"Rank "<<rank<<": assembler memory usage: "<<count<<" KiB"<<endl;
			cout.flush();
			break;
		}
	}
	f.close();
	#endif
}

uint64_t getPathUniqueId(int rank,int id){
	uint64_t a=id;
	a=a*MAX_NUMBER_OF_MPI_PROCESSES+rank;
	return a;
}

int getIdFromPathUniqueId(uint64_t a){
	return a/MAX_NUMBER_OF_MPI_PROCESSES;
}

int getRankFromPathUniqueId(uint64_t a){
	int rank=a%MAX_NUMBER_OF_MPI_PROCESSES;
	return rank;
}

void now(){
	#ifdef OS_POSIX
	time_t m_endingTime=time(NULL);
	struct tm * timeinfo;
	timeinfo=localtime(&m_endingTime);
	cout<<"Date: "<<asctime(timeinfo);
	#endif
}

void print64(uint64_t a){
	for(int k=63;k>=0;k-=2){
		int bit=a<<(k-1)>>63;
		printf("%i",bit);
		bit=a<<(k)>>63;
		printf("%i ",bit);
	}
	printf("\n");
}

void print8(uint8_t a){
	for(int i=7;i>=0;i--){
		int bit=((((uint64_t)a)<<((sizeof(uint64_t)*8-1)-i))>>(sizeof(uint64_t)*8-1));
		printf("%i ",bit);
	}
	printf("\n");
}

uint8_t charToCode(char a){
	switch (a){
		case 'A':
		case '0':
			return RAY_NUCLEOTIDE_A;
		case 'C':
		case '1':
			return RAY_NUCLEOTIDE_C;
		case 'G':
		case '2':
			return RAY_NUCLEOTIDE_G;
		case 'T':
		case '3':
			return RAY_NUCLEOTIDE_T;
		default:
			return RAY_NUCLEOTIDE_A;
	}
}

char codeToChar(uint8_t a,bool color){
	if(color){
		switch(a){
			case RAY_NUCLEOTIDE_A:
				return DOUBLE_ENCODING_A_COLOR;
			case RAY_NUCLEOTIDE_T:
				return DOUBLE_ENCODING_T_COLOR;
			case RAY_NUCLEOTIDE_C:
				return DOUBLE_ENCODING_C_COLOR;
			case RAY_NUCLEOTIDE_G:
				return DOUBLE_ENCODING_G_COLOR;
		}
		return DOUBLE_ENCODING_A_COLOR;
	}

	switch(a){
		case RAY_NUCLEOTIDE_A:
			return 'A';
		case RAY_NUCLEOTIDE_T:
			return 'T';
		case RAY_NUCLEOTIDE_C:
			return 'C';
		case RAY_NUCLEOTIDE_G:
			return 'G';
	}
	return 'A';
}

int portableProcessId(){
	#ifdef OS_POSIX
	return getpid();
	#elif defined(OS_WIN)
	return GetCurrentProcessId();
	#else
	return -1;
	#endif
}
