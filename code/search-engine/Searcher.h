/*
 	Ray
    Copyright (C) 2010, 2011, 2012 Sébastien Boisvert

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

#ifndef _Searcher_h
#define _Searcher_h

#include <core/Parameters.h>
#include <scheduling/SwitchMan.h>
#include <communication/VirtualCommunicator.h>
#include <memory/RingAllocator.h>
#include <search-engine/SearchDirectory.h>
#include <structures/StaticVector.h>
#include <assembler/TimePrinter.h>
#include <stdint.h>
#include <fstream>
#include <stdio.h> /*for FILE */
#include <vector>
#include <string>
#include <communication/BufferedData.h>
#include <profiling/Derivative.h>
#include <map>
using namespace std;

/**
 * This class searches for sequences in the de Bruijn graph
 * It outputs biological abundance readouts
 * \author Sébastien Boisvert
 **/
class Searcher{
	/** number of pending messages */
	int m_pendingMessages;
	
	int m_numberOfRanksThatFinishedSequenceAbundances;

	/** total number of kmers processed */
	uint64_t m_kmersProcessed;

	/** manually buffered data */
	BufferedData m_bufferedData;

	/** point-wise derivative for speed readouts */
	Derivative m_derivative;

	/** the Ray parameters */
	Parameters*m_parameters;

	/** the switchman */
	SwitchMan*m_switchMan;
	VirtualCommunicator*m_virtualCommunicator;
	StaticVector*m_outbox;
	StaticVector*m_inbox;
	RingAllocator*m_outboxAllocator;
	TimePrinter*m_timePrinter;

	/** the identifier to give to the virtual communicator */
	uint64_t m_workerId;

	// for counting entries in category files
	bool m_countElementsMasterStarted;
	bool m_countElementsSlaveStarted;

	// for counting stuff in contigs
	bool m_countContigKmersSlaveStarted;
	bool m_countContigKmersMasterStarted;

	// for counting sequences
	bool m_listedDirectories;
	bool m_countedDirectories;
	bool m_shareCounts;
	bool m_sharedCounts;
	int m_directoryIterator;
	int m_fileIterator;
	bool m_waiting;

	/** option that indicates if detailed reports are needed */
	bool m_writeDetailedFiles;

	/** contig paths */
	vector<vector<Kmer> >*m_contigs;
	vector<uint64_t>*m_contigNames;

	// synchronization
	int m_ranksDoneCounting;
	int m_ranksDoneSharing;
	int m_ranksSynced;
	bool m_synchronizationIsDone;

	int m_masterDirectoryIterator;
	int m_masterFileIterator;
	bool m_sendCounts;

	/** contig iterator */
	int m_contig;

	/** contig position iterator */
	int m_contigPosition;

	/** only used by master, tells which rank to call next */
	int m_rankToCall;

	/** indicates if the rank has sent a control message already */
	bool m_finished;

	// for the virtual communicator
	vector<uint64_t> m_activeWorkers;

	/** file to write coverage distribution */
	ofstream m_currentCoverageFile;
	bool m_waitingForAbundanceReply;

	map<int,int> m_coverageDistribution;

	/** search directory objects */
	SearchDirectory*m_searchDirectories;
	int m_searchDirectories_size;

	/** for the master rank */
	ofstream m_contigSummaryFile;
	ofstream m_identificationFile;

	// base names of directories
	vector<string> m_directoryNames;

	// base names of files
	vector<vector<string> > m_fileNames;

	/** number of matches for a sequence */
	int m_matches;

	/** k-mer length */
	int m_kmerLength;

	/** flag */
	bool m_requestedCoverage;

	int m_numberOfKmers;

	/** files to write */
	map<int,map<int,FILE*> > m_arrayOfFiles;

	/** track the descriptors */
	int m_activeFiles;

	/** partition */
	int m_firstSequence;
	int m_lastSequence;

	/** states for sequence abundances */
	bool m_countSequenceKmersSlaveStarted;
	bool m_countSequenceKmersMasterStarted;

	/** state */
	bool m_createdSequenceReader;

	/** iteraetor */
	int m_sequenceIterator;
	int m_globalSequenceIterator;
	int m_globalFileIterator;


//_-----------------_//


	void showContigAbundanceProgress();
	void createTrees();

	void showSequenceAbundanceProgress();

	bool isSequenceOwner();

	void printDirectoryStart();

	string getBaseName(string a);

	void showProcessedKmers();
public:

	void countElements_masterMethod();
	void countElements_slaveMethod();

	void countContigKmers_masterHandler();
	void countContigKmers_slaveHandler();

	void countSequenceKmers_masterHandler();
	void countSequenceKmers_slaveHandler();
	
	void constructor(Parameters*parameters,StaticVector*outbox,TimePrinter*timePrinter,SwitchMan*switchMan,
	VirtualCommunicator*m_vc,StaticVector*inbox,RingAllocator*outboxAllocator);

	void setContigs(vector<vector<Kmer> >*paths,vector<uint64_t>*names);

}; /* Searcher */

#endif

