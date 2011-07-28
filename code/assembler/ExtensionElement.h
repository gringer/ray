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

#ifndef _ExtensionElement
#define _ExtensionElement

#include <structures/PairedRead.h>
#include <memory/MyAllocator.h>

/*
 * An extension element is a read mapped on an
 * extension for an MPI rank.
 */
class ExtensionElement{
	PairedRead m_pairedRead;
	char*m_readSequence;
	int m_position;
	uint16_t m_strandPosition;
	char m_strand;
	bool m_hasPairedRead;
	uint8_t m_type;
	bool m_canMove;
public:
	void setStrandPosition(int a);
	int getStrandPosition();
	void setStartingPosition(int a);
	void setStrand(char a);
	void setPairedRead(PairedRead a);
	void setType(int a);
	bool hasPairedRead();
	void setSequence(string b,MyAllocator* allocator);
	int getPosition();
	char getStrand();
	PairedRead*getPairedRead();
	char*getSequence();
	bool isLeftEnd();
	bool isRightEnd();
	int getType();
	void removeSequence();
	void constructor();
	bool canMove();
	void freezePlacement();
} ATTRIBUTE_PACKED;

#endif
