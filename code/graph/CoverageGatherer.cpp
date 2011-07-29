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

#include <graph/CoverageGatherer.h>
#ifdef ASSERT
#include <assert.h>
#endif
#include <communication/Message.h>
#include <communication/mpi_tags.h>
#include <core/slave_modes.h>
#include <stdio.h>
#include <stdint.h>
#include <core/constants.h>
#include <structures/Kmer.h>
#include <graph/GridTableIterator.h>
#include <graph/KmerAcademyIterator.h>
#include <sstream>
using namespace std;

void CoverageGatherer::writeKmers(){
	#ifdef ASSERT
	uint64_t n=0;
	#endif
	int m_wordSize = m_parameters->getWordSize();
	if(m_subgraph->getKmerAcademy()->size()==0){
		(*m_slaveMode)=RAY_SLAVE_MODE_DO_NOTHING;
		Message aMessage(NULL,0,MASTER_RANK,RAY_MPI_TAG_COVERAGE_END,
			m_parameters->getRank());
		(*m_outbox).push_back(aMessage);
		return;
	}
	GridTableIterator iterator;
	iterator.constructor(m_subgraph,m_wordSize,m_parameters);
	FILE*kmerFile=NULL;
	ostringstream name;
	name<<m_parameters->getPrefix()<<".kmers.txt";
	if(m_parameters->getRank()==0)
		kmerFile=fopen(name.str().c_str(),"w"); // create empty file 
	else
		kmerFile=fopen(name.str().c_str(),"a"); // append to file 

	fprintf(kmerFile,"# Generated by Ray %s\n",RAY_VERSION);
	fprintf(kmerFile,"# The length of k-mers is %i\n",m_parameters->getWordSize());
	fprintf(kmerFile,"# This file contains the k-mer graph (subgraph of the de Bruijn graph for alphabet={A,C,G,T})\n");
	fprintf(kmerFile,"# The genome sequence you are looking for is a path in this maze\n");
	fprintf(kmerFile,"# You can kickstart your assembly algorithm development by loading this file\n");
	fprintf(kmerFile,"# Note that vertices with a coverage of 1 are not considered.\n");
	fprintf(kmerFile,"# Format:\n");
	fprintf(kmerFile,"# k-mer sequence; coverage value; first nucleotide of parents; last nucleotide of children\n");
	fprintf(kmerFile,"# Example in nucleotide space:\n# ATCG;10;T G;C A\n#  ATCG has a coverage value of 10\n");
	fprintf(kmerFile,"#  Ingoing arcs: TATC -> ATCG and GATC -> ATCG\n");
	fprintf(kmerFile,"#  Outgoing arcs: ATCG -> TCGC and ATCG -> TCGA\n");

	while(iterator.hasNext()){
		Vertex*node=iterator.next();
		Kmer key=*(iterator.getKey());
		int coverage=node->getCoverage(&key);
		m_distributionOfCoverage[coverage]++;
		#ifdef ASSERT
		n++;
		#endif
		string kmerSequence=key.toBSString(m_wordSize);
		vector<Kmer> parents=node->getIngoingEdges(&key,m_wordSize);
		vector<Kmer> children=node->getOutgoingEdges(&key,m_wordSize);
		fprintf(kmerFile,"%s;%i;",kmerSequence.c_str(),coverage);
		for(int i=0;i<(int)parents.size();i++){
			char edge=parents[i].getFirstSymbol(m_wordSize, false);
			if(i!=0)
				fprintf(kmerFile," ");

			fprintf(kmerFile,"%c",edge);
		}
		fprintf(kmerFile,";");
		for(int i=0;i<(int)children.size();i++){
			char edge=children[i].getLastSymbol(m_wordSize, false);
			if(i!=0)
				fprintf(kmerFile," ");

			fprintf(kmerFile,"%c",edge);
		}
		fprintf(kmerFile,"\n");
	}
	fclose(kmerFile);
		
	#ifdef ASSERT
	if(n!=m_subgraph->size()){
		cout<<"n="<<n<<" size="<<m_subgraph->size()<<endl;
	}
	assert(n==m_subgraph->size());
	#endif
	m_waiting=false;
	m_coverageIterator=m_distributionOfCoverage.begin();
}

void CoverageGatherer::work(){

	if(m_distributionOfCoverage.size()==0){
		#ifdef ASSERT
		uint64_t n=0;
		#endif
		if(m_subgraph->getKmerAcademy()->size()==0){
			(*m_slaveMode)=RAY_SLAVE_MODE_DO_NOTHING;
			Message aMessage(NULL,0,MASTER_RANK,RAY_MPI_TAG_COVERAGE_END,
				m_parameters->getRank());
			(*m_outbox).push_back(aMessage);
			return;
		}
		KmerAcademyIterator iterator;
		iterator.constructor(m_subgraph->getKmerAcademy(),m_parameters->getWordSize(),m_parameters);
		while(iterator.hasNext()){
			KmerCandidate*node=iterator.next();
			Kmer key=*(iterator.getKey());
			int coverage=node->m_count;
			m_distributionOfCoverage[coverage]++;
			#ifdef ASSERT
			n++;
			#endif
		}
			
		#ifdef ASSERT
		if(n!=m_subgraph->getKmerAcademy()->size()){
			cout<<"Expected (from iterator)="<<n<<" Actual (->size())="<<m_subgraph->getKmerAcademy()->size()<<endl;
		}
		assert(n==m_subgraph->getKmerAcademy()->size());
		#endif
		m_waiting=false;
		m_coverageIterator=m_distributionOfCoverage.begin();
	}else if(m_waiting){
		if((*m_inbox).size()>0&&(*m_inbox)[0]->getTag()==RAY_MPI_TAG_COVERAGE_DATA_REPLY){
			m_waiting=false;
		}
	}else{
		uint64_t*messageContent=(uint64_t*)m_outboxAllocator->allocate(MAXIMUM_MESSAGE_SIZE_IN_BYTES);
		int count=0;
		int maximumElements=MAXIMUM_MESSAGE_SIZE_IN_BYTES/sizeof(uint64_t);
		while(count<maximumElements && m_coverageIterator!=m_distributionOfCoverage.end()){
			int coverage=m_coverageIterator->first;
			uint64_t numberOfVertices=m_coverageIterator->second;
			messageContent[count]=coverage;
			messageContent[count+1]=numberOfVertices;
			count+=2;
			m_coverageIterator++;
		}

		if(count!=0){
			Message aMessage(messageContent,count,MASTER_RANK,RAY_MPI_TAG_COVERAGE_DATA,
				m_parameters->getRank());
			
			(*m_outbox).push_back(aMessage);
			m_waiting=true;
		}else{
			m_distributionOfCoverage.clear();
			(*m_slaveMode)=RAY_SLAVE_MODE_DO_NOTHING;
			Message aMessage(NULL,0,MASTER_RANK,RAY_MPI_TAG_COVERAGE_END,
				m_parameters->getRank());
			(*m_outbox).push_back(aMessage);
		}
	}
}

void CoverageGatherer::constructor(Parameters*parameters,StaticVector*inbox,StaticVector*outbox,int*slaveMode,
	GridTable*subgraph,RingAllocator*outboxAllocator){
	m_parameters=parameters;
	m_slaveMode=slaveMode;
	m_outboxAllocator=outboxAllocator;
	m_inbox=inbox;
	m_outbox=outbox;
	m_subgraph=subgraph;
}

