/*
 	Ray
    Copyright (C) 2011, 2012  Sébastien Boisvert

	http://DeNovoAssembler.SourceForge.Net/

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You have received a copy of the GNU General Public License
    along with this program (gpl-3.0.txt).  
	see <http://www.gnu.org/licenses/>

*/

#ifndef _KmerAcademy_H
#define _KmerAcademy_H

#include <memory/MyAllocator.h>
#include <plugin_KmerAcademyBuilder/Kmer.h>
#include <plugin_KmerAcademyBuilder/KmerCandidate.h>
#include <application_core/Parameters.h>
#include <structures/MyHashTable.h>
#include <stdint.h>


/**
 * The KmerAcademy is the place where KmerCandidate  
 * train to become  Vertex.
 * If they fail to endure their hard training,
 * they remain here. Usually, those failing are observed
 * only once.
 * The KmerAcademy is burned and destroyed once the graph is ready.
 * \author Sébastien Boisvert
 */
class KmerAcademy{
	Parameters*m_parameters;
	LargeCount m_size;
	bool m_inserted;
	MyHashTable<Kmer,KmerCandidate> m_hashTable;

public:
	void constructor(Rank rank,Parameters*a);
	LargeCount size();
	KmerCandidate*find(Kmer*key);
	KmerCandidate*insert(Kmer*key);
	bool inserted();
	void destructor();
	MyHashTable<Kmer,KmerCandidate>*getHashTable();
	void printStatistics();
	void completeResizing();
};

#endif
