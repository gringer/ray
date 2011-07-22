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

#include <core/constants.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assembler/KmerAcademyBuilder.h>
#include <assert.h>
#include <communication/Message.h>
#include <time.h>
#include <structures/StaticVector.h>
#include <core/common_functions.h>
#include <memory/malloc_types.h>

void KmerAcademyBuilder::process(int*m_mode_send_vertices_sequence_id,
				ArrayOfReads*m_myReads,
				bool*m_reverseComplementVertex,
				int rank,
				StaticVector*m_outbox,
				StaticVector*m_inbox,
				bool*m_mode_send_vertices,
				int wordSize,
				int size,
				RingAllocator*m_outboxAllocator,
				int*m_mode
				){
	if(this->m_outbox==NULL){
		m_rank=rank;
		this->m_mode=m_mode;
		this->m_outbox=m_outbox;
		this->m_outboxAllocator=m_outboxAllocator;
	}
	#ifdef ASSERT
	assert(m_pendingMessages>=0);
	#endif
	if(m_inbox->size()>0&&m_inbox->at(0)->getTag()==RAY_MPI_TAG_KMER_ACADEMY_DATA_REPLY){
		m_pendingMessages--;
	}

	if(m_pendingMessages!=0){
		return;
	}

	if(m_finished)
		return;

	if(*m_mode_send_vertices_sequence_id%10000==0 &&m_mode_send_vertices_sequence_id_position==0
	&&*m_mode_send_vertices_sequence_id<(int)m_myReads->size()){
		string reverse="";
		if(*m_reverseComplementVertex==true){
			reverse="(reverse complement) ";
		}
		printf("Rank %i is counting k-mers in sequence reads %s[%i/%i]\n",rank,reverse.c_str(),(int)*m_mode_send_vertices_sequence_id+1,(int)m_myReads->size());
		fflush(stdout);
	}

	if(*m_mode_send_vertices_sequence_id>(int)m_myReads->size()-1){
		// flush data
		flushAll(m_outboxAllocator,m_outbox,rank);
		if(m_pendingMessages==0){
			#ifdef ASSERT
			assert(m_bufferedData.isEmpty());
			#endif

			Message aMessage(NULL,0, MASTER_RANK, 
				RAY_MPI_TAG_KMER_ACADEMY_DISTRIBUTED,rank);
			m_outbox->push_back(aMessage);
			m_finished=true;
			printf("Rank %i is counting k-mers in sequence reads [%i/%i] (completed)\n",rank,(int)*m_mode_send_vertices_sequence_id,(int)m_myReads->size());
			fflush(stdout);
			m_bufferedData.showStatistics(m_parameters->getRank());
		}
	}else{
		if(m_mode_send_vertices_sequence_id_position==0){
			m_readSequence = (*m_myReads)[(*m_mode_send_vertices_sequence_id)]->getSeq(true,false);
		
		}
		int len=m_readSequence.length();

		if(len<wordSize){
			(*m_mode_send_vertices_sequence_id)++;
			(m_mode_send_vertices_sequence_id_position)=0;
			return;
		}

		int lll=len-wordSize+1;
		
		#ifdef ASSERT
		assert(m_readSequence.length()>0);
		#endif

		int p=(m_mode_send_vertices_sequence_id_position);
		Kmer a(m_readSequence,p,wordSize);
		if(a.isValid()){

			// reverse complement
			Kmer aRC=a.rComp(wordSize);

			bool aLower = a.isLower(aRC);

			int rankToFlush=(aLower?a:aRC).hash_function_1()%m_parameters->getSize();
			
			for(int i=0;i<KMER_U64_ARRAY_SIZE;i++){
				m_bufferedData.addAt(rankToFlush,(aLower?a:aRC).getRawBits(i));
			}

			if(m_bufferedData.flush(rankToFlush,KMER_U64_ARRAY_SIZE,RAY_MPI_TAG_KMER_ACADEMY_DATA,m_outboxAllocator,m_outbox,rank,false)){
				m_pendingMessages++;
			}
		}

		(m_mode_send_vertices_sequence_id_position++);
		if((m_mode_send_vertices_sequence_id_position)==lll){
			(*m_mode_send_vertices_sequence_id)++;
			(m_mode_send_vertices_sequence_id_position)=0;
		}
	}
}

void KmerAcademyBuilder::constructor(int size,Parameters*parameters,GridTable*graph){
	m_subgraph=graph;
	m_parameters=parameters;
	m_finished=false;
	m_distributionIsCompleted=false;
	m_outbox=NULL;
	m_mode_send_vertices_sequence_id_position=0;
	m_bufferedData.constructor(size,MAXIMUM_MESSAGE_SIZE_IN_BYTES/sizeof(uint64_t),RAY_MALLOC_TYPE_KMER_ACADEMY_BUFFER,m_parameters->showMemoryAllocations(),KMER_U64_ARRAY_SIZE);
	
	m_pendingMessages=0;
	m_size=size;
}

void KmerAcademyBuilder::setReadiness(){
	#ifdef ASSERT
	assert(m_pendingMessages>0);
	#endif
	m_pendingMessages--;
}

void KmerAcademyBuilder::flushAll(RingAllocator*m_outboxAllocator,StaticVector*m_outbox,int rank){
	if(!m_bufferedData.isEmpty()){
		m_pendingMessages+=m_bufferedData.flushAll(RAY_MPI_TAG_KMER_ACADEMY_DATA,m_outboxAllocator,m_outbox,rank);
		return;
	}
}

bool KmerAcademyBuilder::finished(){
	return m_finished;
}

void KmerAcademyBuilder::assertBuffersAreEmpty(){
	assert(m_bufferedData.isEmpty());
	assert(m_mode_send_vertices_sequence_id_position==0);
}

void KmerAcademyBuilder::incrementPendingMessages(){
	m_pendingMessages++;
}

void KmerAcademyBuilder::showBuffers(){
	#ifdef ASSERT
	assert(m_bufferedData.isEmpty());
	#endif
}

bool KmerAcademyBuilder::isDistributionCompleted(){
	return m_distributionIsCompleted;
}

void KmerAcademyBuilder::setDistributionAsCompleted(){
	m_distributionIsCompleted=true;
}
