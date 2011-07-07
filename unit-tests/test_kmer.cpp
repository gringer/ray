
#include <structures/Vertex.h>
#include <unit-tests/unitTest.h>
#include <structures/Kmer.h>
#include <core/constants.h>
#include <set>
#include <core/common_functions.h>
#include <vector>
#include <assert.h>
#include <string>
#include <iostream>
using namespace std;

void test_addInEdge(){
	string a="AGCAAGTTAGCAACATCATATGAGTGCAATCCTGTTGTAGGCTCATCTAAGACATAAATAGTT";
	string b= "GCAAGTTAGCAACATCATATGAGTGCAATCCTGTTGTAGGCTCATCTAAGACATAAATAGTTT";
	int wordSize=a.length();

	Kmer aKmer(a);
	Kmer bKmer(b);

	Vertex bVertex;
	bVertex.constructor();
	bVertex.m_lowerKey=bKmer;
	bVertex.addIngoingEdge(&bKmer,&aKmer,wordSize);
	
	vector<Kmer>inEdges=bVertex.getIngoingEdges(&bKmer,wordSize);
	bool found=false;
	for(int j=0;j<(int)inEdges.size();j++){
		if(inEdges[j]==aKmer){
			found=true;
			break;
		}
	}
	if(!found){
		cout<<"Expected: "<<a<<endl;
		cout<<"Actual:"<<endl;
		cout<<inEdges.size()<<endl;
		for(int j=0;j<(int)inEdges.size();j++){
			cout<<inEdges[j].toString(wordSize, true) << " ("
					<< inEdges[j].toBSString(wordSize) << ")"<<endl;
		}
	}
	assertEquals(inEdges.size(),1);
	assertEquals(found,true);
}

void test_addOutEdge(){
	string a="CAATAAGTAAAAAAGATTTTGTAACTTTCACAGCCTTATTTTTATCAATAGATACTGATAT";
	string b= "AATAAGTAAAAAAGATTTTGTAACTTTCACAGCCTTATTTTTATCAATAGATACTGATATT";
	int wordSize=a.length();

	Kmer aKmer(a);
	Kmer bKmer(b);

	Vertex aVertex;
	aVertex.constructor();
	Kmer lower=aKmer;
	Kmer aRC=aKmer.rComp(wordSize);

	if(aRC<lower){
		lower=aRC;
	}
	aVertex.m_lowerKey=lower;
	aVertex.addOutgoingEdge(&aKmer,&bKmer,wordSize);
	
	vector<Kmer>Edges=aVertex.getOutgoingEdges(&aKmer,wordSize);
	bool found=false;
	for(int j=0;j<(int)Edges.size();j++){
		if(Edges[j]==bKmer){
			found=true;
			break;
		}
	}
	if(!found){
		cout<<"Expected: "<<endl;
		cout<<b<<endl;
		cout<<"Actual:"<<endl;
		for(int j=0;j<(int)Edges.size();j++){
			cout<<Edges[j].toString(wordSize, true) << " ("
					<< Edges[j].toBSString(wordSize) << ")"<<endl;
		}
		uint8_t edges=aVertex.getEdges(&aKmer);
		cout<<"Edges"<<endl;
		print8(edges);
	}
	assertEquals(Edges.size(),1);
	assertEquals(found,true);
}

void test_addInEdge2(){
	string a="AGCAAGTTAGCAACATCATATGAGTGCAATCCTGTTGTAGGCTCATCTAAGACATAAATAGTT";
	string b= "GCAAGTTAGCAACATCATATGAGTGCAATCCTGTTGTAGGCTCATCTAAGACATAAATAGTTT";
	int wordSize=a.length();

	Kmer aKmer(a);
	Kmer bKmer(b);

	Vertex bVertex;
	bVertex.constructor();
	Kmer bRC=bKmer.rComp(wordSize);
	Kmer lower=bKmer;
	Kmer aRC=aKmer.rComp(wordSize);

	if(bRC<lower){
		lower=bRC;
	}
	bVertex.m_lowerKey=lower;
	bVertex.addIngoingEdge(&bKmer,&aKmer,wordSize);
	
	vector<Kmer>inEdges=bVertex.getIngoingEdges(&bKmer,wordSize);
	bool found=false;
	for(int j=0;j<(int)inEdges.size();j++){
		if(inEdges[j]==aKmer){
			found=true;
			break;
		}
	}
	if(!found){
		cout<<"Expected: "<<a<<endl;
		cout<<"Actual:"<<endl;
		cout<<inEdges.size()<<endl;
		for(int j=0;j<(int)inEdges.size();j++){
			cout<<inEdges[j].toString(wordSize, true) << " ("
					<< inEdges[j].toBSString(wordSize) << ")"<<endl;
		}
	}
	assertEquals(inEdges.size(),1);
	assertEquals(found,true);
}

void test_out_large(){
	string a="TCAAAAATTTCTTTCAAAGTAATCTCATAAGCTGCTGGA";
	string b= "CAAAAATTTCTTTCAAAGTAATCTCATAAGCTGCTGGAT";
	int wordSize=a.length();

	// description of m_edges:
	// outgoing  ingoing
	//
	// G C T A G C T A
	//
	// 7 6 5 4 3 2 1 0
	
	uint8_t edges=(1<<(4+RAY_NUCLEOTIDE_T));

	Kmer aKmer(a);
	Kmer bKmer(b);
	
	vector<Kmer>oEdges=aKmer.getOutgoingEdges(edges,wordSize);
	assertEquals(oEdges.size(),1);

	Kmer actual=oEdges[0];
	string actualStr=actual.toBSString(wordSize);
	if(actualStr!=b){
		cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
		cout<<"WordSize: "<<wordSize<<endl;
		cout<<"Expected"<<endl;
		cout<<a<<" -> "<<b<<endl;
		cout<<"Actual:"<<endl;
		cout<<a<<" -> "<<actualStr<<"*"<<endl;
		cout<<endl;
	}
	assertEquals(actualStr,b);

	if(actual!=bKmer){
		cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
		cout<<"WordSize: "<<wordSize<<endl;
		cout<<"Expected: "<<endl;
		bKmer.printPieces();
		cout<<"Actual: "<<endl;
		actual.printPieces();
	}
	assert(actual==bKmer);
}

void test_Ingoing_large2(){
	string a="AGCAAGTTAGCAACATCATATGAGTGCAATCCTGTTGTAGGCTCATCTAAGACATAAATAGTT";
	string b= "GCAAGTTAGCAACATCATATGAGTGCAATCCTGTTGTAGGCTCATCTAAGACATAAATAGTTT";
	int wordSize=a.length();

	// description of m_edges:
	// outgoing  ingoing
	//
	// G C T A G C T A
	//
	// 7 6 5 4 3 2 1 0
	
	uint8_t edges=(1<<0);

	Kmer aKmer(a);
	assertEquals(aKmer.getFirstCode(false),RAY_NUCLEOTIDE_A);

	Kmer bKmer(b);
	vector<Kmer>inEdges=bKmer.getIngoingEdges(edges,wordSize);
	Kmer actual=inEdges[0];
	string actualStr=actual.toBSString(wordSize);
	if(actualStr!=a){
		cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
		cout<<"WordSize: "<<wordSize<<endl;
		cout<<"Expected"<<endl;
		cout<<a<<" -> "<<b<<endl;
		cout<<"Actual:"<<endl;
		cout<<actualStr<<" -> "<<b<<endl;
		cout<<endl;
	}
	assertEquals(actualStr,a);

	if(actual!=aKmer){
		cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
		cout<<"WordSize: "<<wordSize<<endl;
		cout<<"Expected: "<<endl;
		aKmer.printPieces();
		cout<<"Actual: "<<endl;
		actual.printPieces();
	}
	assert(actual==aKmer);
}




void test_Ingoing_large(){
	string a="TCAAAAATTTCTTTCAAAGTAATCTCATAAGCTGCTGGA";
	string b= "CAAAAATTTCTTTCAAAGTAATCTCATAAGCTGCTGGAT";
	int wordSize=a.length();

	// description of m_edges:
	// outgoing  ingoing
	//
	// G C T A G C T A
	//
	// 7 6 5 4 3 2 1 0
	
	uint8_t edges=(1<<RAY_NUCLEOTIDE_T);

	Kmer aKmer(a);
	Kmer bKmer(b);
	
	vector<Kmer>inEdges=bKmer.getIngoingEdges(edges,wordSize);
	Kmer actual=inEdges[0];
	string actualStr=actual.toBSString(wordSize);
	if(actualStr!=a){
		cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
		cout<<"WordSize: "<<wordSize<<endl;
		cout<<"Expected"<<endl;
		cout<<a<<" -> "<<b<<endl;
		cout<<"Actual:"<<endl;
		cout<<actualStr<<" -> "<<b<<endl;
		cout<<endl;
	}
	assertEquals(actualStr,a);

	if(actual!=aKmer){
		cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
		cout<<"WordSize: "<<wordSize<<endl;
		cout<<"Expected: "<<endl;
		aKmer.printPieces();
		cout<<"Actual: "<<endl;
		actual.printPieces();
	}
	assert(actual==aKmer);
}



void test_out(){
	string a="GACTTGATTAGACAAGAAGTT";
	string b= "ACTTGATTAGACAAGAAGTTG";
	int wordSize=a.length();

	// description of m_edges:
	// outgoing  ingoing
	//
	
	uint8_t edges=(1<<(4+RAY_NUCLEOTIDE_G));

	Kmer aKmer(a);
	Kmer bKmer(b);
	
	vector<Kmer>oEdges=aKmer.getOutgoingEdges(edges,wordSize);
	assertEquals(oEdges.size(),1);

	Kmer actual=oEdges[0];
	string actualStr=actual.toBSString(wordSize);
	if(actualStr!=b){
		cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
		cout<<"WordSize: "<<wordSize<<endl;
		cout<<"Expected"<<endl;
		cout<<a<<" -> "<<b<<endl;
		cout<<"Actual:"<<endl;
		cout<<a<<" -> "<<actualStr<<"*"<<endl;
		cout<<endl;
	}
	assertEquals(actualStr,b);

	if(actual!=bKmer){
		cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
		cout<<"WordSize: "<<wordSize<<endl;
		cout<<"Expected: "<<endl;
		bKmer.printPieces();
		cout<<"Actual: "<<endl;
		actual.printPieces();
	}
	assert(actual==bKmer);
}

void test_Ingoing(){
	string a="GACTTGATTAGACAAGAAGTT";
	string b= "ACTTGATTAGACAAGAAGTTG";
	int wordSize=a.length();

	// description of m_edges:
	// outgoing  ingoing
	//
	// T G C A T G C A
	//
	// 7 6 5 4 3 2 1 0
	
	uint8_t edges=(1<<RAY_NUCLEOTIDE_G);

	Kmer aKmer(a);
	Kmer bKmer(b);
	
	vector<Kmer>inEdges=bKmer.getIngoingEdges(edges,wordSize);
	Kmer actual=inEdges[0];
	string actualStr=actual.toBSString(wordSize);
	if(actualStr!=a){
		cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
		cout<<"WordSize: "<<wordSize<<endl;
		cout<<"Expected"<<endl;
		cout<<a<<" -> "<<b<<endl;
		cout<<"Actual:"<<endl;
		cout<<actualStr<<" -> "<<b<<endl;
		cout<<endl;
	}
	assertEquals(actualStr,a);

	if(actual!=aKmer){
		cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
		cout<<"WordSize: "<<wordSize<<endl;
		cout<<"Expected: "<<endl;
		aKmer.printPieces();
		cout<<"Actual: "<<endl;
		actual.printPieces();
	}
	assert(actual==aKmer);
}


int main(int argc,char**argv){
	string seq=argv[1];
	int wordSize=seq.length();

	Kmer id(seq);
	Kmer empty;
	string result=id.toBSString(wordSize);
	assert(seq==result);
	char last=seq.at(seq.length()-1);
	char observed=id.getLastSymbol(wordSize, false);
	assert(observed==last);


	// reverse complement
	string rc=reverseComplement(&seq);

	Kmer comp=id.rComp(wordSize);

	string result2=comp.toBSString(wordSize);
	assertEquals(rc,result2);

	Kmer rcId(rc);
	assert(rcId==comp);
	
	// ingoing edges.
	//
	uint8_t edges=0xff;
	vector<Kmer>inEdges=id.getIngoingEdges(edges,wordSize);
	assertEquals(4,inEdges.size());
	set<string> tmp;
	for(int i=0;i<(int)inEdges.size();i++){
		Kmer theKmer=inEdges[i];
		string a=theKmer.toBSString(wordSize);
		Kmer id(a);

		if(theKmer!=id){
			cout<<"MAXKMERLENGTH: "<<MAXKMERLENGTH<<endl;
			cout<<"WordSize: "<<wordSize<<endl;
			cout<<"Expected: "<<endl;
			id.printPieces();
			cout<<"Actual: "<<endl;
			theKmer.printPieces();
		}

		assert(id==theKmer);

		assertEquals(tmp.count(a),0);
		tmp.insert(a);
		assertEquals(a.substr(1,wordSize-1),seq.substr(0,wordSize-1));
	}
	assertEquals(tmp.size(),4);

	// test outgoing edges
	vector<Kmer> outEdges=id.getOutgoingEdges(edges,wordSize);
	assertEquals(4,outEdges.size());
	tmp.clear();
	for(int i=0;i<(int)outEdges.size();i++){
		Kmer theKmer=outEdges[i];
		string a=theKmer.toBSString(wordSize);

		Kmer id(a);// make sure that all bit are set to 0 except those relevant
		if(theKmer!=id){
			cout<<"Expected: "<<endl;
			id.printPieces();
			cout<<"Actual: "<<endl;
			theKmer.printPieces();
		}
		assert(theKmer==id);
		assertEquals(tmp.count(a),0);
		tmp.insert(a);
		assertEquals(seq.substr(1,wordSize-1),a.substr(0,wordSize-1));
	}

	assertEquals(tmp.size(),4);

	test_Ingoing();
	test_out();

	if(MAXKMERLENGTH>32){
		test_Ingoing_large();
		test_Ingoing_large2();
		test_out_large();
		test_addInEdge();
		test_addInEdge2();
		test_addOutEdge();
	}
	return 0;
}

