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

#ifndef _KmerCandidate_H
#define _KmerCandidate_H

#include "plugin_KmerAcademyBuilder/Kmer.h"

/** a KmerCandidate may become a genomic super-star --
 * one of the fews that make it to the graph
 * \author Sébastien Boisvert
 */
class KmerCandidate{
public:
	Kmer m_lowerKey;
	CoverageDepth m_count;

	Kmer getKey();
	void setKey(Kmer key);

}ATTRIBUTE_PACKED;


#endif /* _KmerCandidate_H */


